"""Published β₂AR inputs and evaluation of the six-nanobody reference."""

import ast
import gzip
import hashlib
import io
import itertools
import json
from importlib.resources import files

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .modeling_backends import extract_appraise_job_name
from .score_calculation import calculate_scores, get_APPRAISE_ranking


def resource(name):
    return files("appraise").joinpath("examples", "b2ar", name)


def panel():
    return json.loads(resource("panel.json").read_text(encoding="utf-8"))


def ordered_inputs(names, self_pairs):
    return {
        f"b2AR_and_{left}_vs_{right}": (left, right)
        for left, right in itertools.product(names, repeat=2)
        if self_pairs or left != right
    }


def reference_measurements():
    content = gzip.decompress(resource("measurements.csv.gz").read_bytes())
    if hashlib.sha256(content).hexdigest() != panel()["measurement_sha256"]:
        raise ValueError("The bundled historical measurements failed their source checksum.")
    return pd.read_csv(io.BytesIO(content))


def evaluate(frame, names, predictions, *, self_pairs, historical=False):
    """Rank a complete panel, preserving both candidate perspectives of each structure."""
    spec = panel()
    jobs = ordered_inputs(names, self_pairs)
    df = frame.copy()
    if df.empty:
        raise ValueError("No nanobody measurements were found. Check the structure directory.")
    if df[["model_name", "peptide_chain", "peptide_name", "peptide_seq", "competitors"]].isna().any().any():
        raise ValueError("Missing model, chain or candidate identity in the measurements.")
    if set(df["receptor_name"]) != {"b2AR"} or set(df["anchor_site"]) != {"N-term"}:
        raise ValueError("These measurements must use the manuscript β₂AR construct and N-terminal anchor.")
    df["job"] = df["model_name"].map(extract_appraise_job_name)
    if set(df["job"]) != set(jobs):
        raise ValueError("Incomplete or unexpected ordered competition inputs.")
    if df.duplicated(["model_name", "peptide_chain"]).any():
        raise ValueError("Duplicate candidate perspective for a model and chain.")
    models = df.groupby("model_name", sort=False)
    if not models.size().eq(2).all() or not df.groupby("job")["model_name"].nunique().eq(predictions).all():
        raise ValueError("Incomplete model coverage. Each ordered input needs all predictions and both perspectives.")
    for model, observations in models:
        expected = sorted(jobs[observations["job"].iloc[0]])
        if sorted(observations["peptide_name"]) != expected:
            raise ValueError(f"Candidate identities disagree with the ordered input for {model}.")
        for row in observations.itertuples():
            seq = spec["candidates"][row.peptide_name]
            # The historical PyMOL writer saved just the first 70-character FASTA line.
            if row.peptide_seq != (seq[:70] if historical else seq):
                raise ValueError(f"Sequence mismatch for {row.peptide_name} in {model}.")
            competitor = ast.literal_eval(row.competitors)
            if not isinstance(competitor, list) or len(competitor) != 1 or sorted([row.peptide_name, *competitor]) != expected:
                raise ValueError(f"Competitor does not match the input for {model}.")
    df["competitor"] = df["competitors"].map(lambda value: ast.literal_eval(value)[0])
    df["receptor_Rminor"] = spec["receptor"]["rminor_angstrom"]
    scored = calculate_scores(df, version=1.2)
    if not np.isfinite(scored["Delta_B"]).all():
        raise ValueError("One or more structures have an undefined APPRAISE score.")
    means = scored.groupby(["peptide_name", "competitor"], observed=True)["Delta_B"].mean().reset_index()
    # Each non-self matchup receives one perspective from each of the reciprocal models.
    off_diagonal = df.loc[df["peptide_name"] != df["competitor"]]
    counts = off_diagonal.groupby(["peptide_name", "competitor"])["model_name"].nunique()
    if not counts.eq(2 * predictions).all():
        raise ValueError("Reciprocal matchups have incomplete distinct-model coverage.")
    ranked = get_APPRAISE_ranking(means, number_of_repeats=2 * predictions)
    summary = ranked.groupby("peptide_name")[["match_points", "match_points_tie_breaker", "Delta_B"]].mean()
    columns = ["match_points", "match_points_tie_breaker", "Delta_B"]
    summary = summary.sort_values(columns, ascending=False, kind="stable").reset_index()
    summary["rank"] = summary.groupby(columns, sort=False)["peptide_name"].transform(
        lambda group: float(np.mean(group.index + 1))
    )
    opponents = ranked.groupby("peptide_name")["competitor"].nunique()
    summary["win_fraction"] = 0.5 + summary["match_points"] / (2 * summary["peptide_name"].map(opponents))
    return scored, ranked.pivot(index="peptide_name", columns="competitor", values="Delta_B"), summary


def comparison(ranking):
    """Apply the predeclared scientific criterion to a complete six-candidate ranking."""
    spec = panel()
    truth = spec["experimental_order"]
    if len(ranking) != len(truth) or set(ranking["peptide_name"]) != set(truth):
        raise ValueError("The scientific comparison requires all six published nanobodies.")
    ranked = ranking.set_index("peptide_name")["rank"].reindex(truth)
    rho = float(spearmanr(range(1, 7), ranked).statistic) if ranked.nunique() > 1 else None
    first = ranked.index[ranked == 1].tolist()
    last = ranked.index[ranked == 6].tolist()
    meets = rho is not None and rho >= np.nextafter(29 / 35, -np.inf) and first == ["Nb6B9"] and last == ["Nb80"]
    return {
        "experimental_order": truth, "historical_order": spec["historical_order"],
        "ranks": {name: float(rank) for name, rank in ranked.items()},
        "spearman_rho": rho, "historical_rho": 31 / 35,
        "delta_from_historical": None if rho is None else rho - 31 / 35,
        "minimum_rho": 29 / 35, "meets_ranking_criterion": bool(meets),
    }
