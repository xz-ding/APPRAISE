import gzip
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from appraise import nanobody_models, runner
from appraise.nanobody import comparison, evaluate, panel, reference_measurements, resource
from appraise.utilities import rank_tournament_results


def read(path):
    return json.loads(path.read_text())


def test_reference_cli_and_python_reproduce_published_result(tmp_path):
    result = runner.run(runner.RunConfig(example="b2ar-reference", output_dir=tmp_path / "api"))
    assert result["status"] == "complete", result
    assert (result["completed_jobs"], result["completed_predictions"], result["completed_measurements"]) == (36, 180, 360)
    assert result["comparison"]["spearman_rho"] == pytest.approx(31 / 35)
    assert result["comparison"]["meets_ranking_criterion"]
    assert pd.read_csv(tmp_path / "api/ranking.csv").peptide_name.tolist() == panel()["historical_order"]
    proc = subprocess.run([sys.executable, "-m", "appraise", "score", "--example", "b2ar-reference",
                           "--output", str(tmp_path / "cli")], capture_output=True, text=True, check=True)
    cli = json.loads(proc.stdout)
    assert result["comparison"] == cli["comparison"]
    assert result["artifacts"]["ranking.csv"]["sha256"] == cli["artifacts"]["ranking.csv"]["sha256"]


def test_missing_duplicate_and_displaced_perspectives():
    source = reference_measurements()
    names = list(panel()["candidates"])
    for data in (source.iloc[1:], pd.concat([source, source.tail(1)]),
                 source.assign(peptide_seq=source.peptide_seq.str[:-1])):
        with pytest.raises(ValueError):
            evaluate(data, names, 5, self_pairs=True, historical=True)


def test_gate_boundary_and_tied_endpoints():
    names = panel()["experimental_order"]
    ranks = pd.DataFrame({"peptide_name": names, "rank": [1, 3, 4, 2, 5, 6]})
    assert comparison(ranks)["spearman_rho"] == pytest.approx(29 / 35)
    assert comparison(ranks)["meets_ranking_criterion"]
    ranks["rank"] = [1.5, 1.5, 3, 4, 5, 6]
    assert not comparison(ranks)["meets_ranking_criterion"]
    ranks["rank"] = [1, 2, 3, 4, 5.5, 5.5]
    assert not comparison(ranks)["meets_ranking_criterion"]
    ranks["rank"] = [3.5] * 6
    assert comparison(ranks)["spearman_rho"] is None


def test_tie_break_interval_includes_both_ends():
    frame = pd.DataFrame({"peptide_name": list("ABC"), "value": [-2., 0., 2.]})
    rank_tournament_results(frame, "value", tie_threshold=1, points=[3, 2, 0], number_of_repeats=2)
    assert frame["match_points"].tolist() == [0, 2, 3]
    assert frame["match_points_tie_breaker"].tolist() == [2, 2, 2]


def test_partial_saved_run_recovers_by_identity(tmp_path, monkeypatch):
    pytest.importorskip("pymol")
    config = runner.RunConfig(example="b2ar-saved", output_dir=str(tmp_path / "saved"), predictions_per_input=1)
    def interrupt(event):
        if event["stage"] == "measured" and event["completed_predictions"] == 1:
            raise KeyboardInterrupt
    result = runner.run(config, progress=interrupt)
    assert result["status"] == "interrupted", result
    root = Path(config.output_dir)
    first_file = next(root.glob("models/*/*/*.pdb"))
    first_mtime = first_file.stat().st_mtime_ns
    first_measure = first_file.with_name("measurements.csv")
    measurement_time = first_measure.stat().st_mtime_ns
    resumed = runner.resume(root)
    assert resumed["status"] == "complete", resumed
    assert (resumed["completed_jobs"], resumed["completed_predictions"], resumed["completed_measurements"]) == (2, 2, 4)
    assert first_file.stat().st_mtime_ns == first_mtime
    assert first_measure.stat().st_mtime_ns == measurement_time
    rows = pd.read_csv(root / "measurements.csv")
    assert rows["peptide_seq"].str.len().eq(120).all()
    assert not rows.duplicated(["model_name", "peptide_chain"]).any()
    _, grid, ranking = evaluate(rows, ["Nb6B9", "Nb80"], 1, self_pairs=False)
    assert np.isnan(np.diag(grid)).all()
    assert ranking["win_fraction"].between(0, 1).all()
    before = (root / "result.json").read_bytes()
    wrong = runner.run(runner.RunConfig(example="b2ar-reference", output_dir=str(root)))
    assert wrong["failure_code"] == "config_mismatch" and (root / "result.json").read_bytes() == before
    first_measure.write_text("corrupt")
    assert runner.resume(root)["status"] == "complete"
    assert first_file.stat().st_mtime_ns == first_mtime
    assert pd.read_csv(first_measure).shape[0] == 2
    manifest_path = root / "run_manifest.json"
    manifest = read(manifest_path)
    manifest["samples"]["b2AR_and_Nb6B9_vs_Nb80/1"]["observed_chains"] = ["C", "C", "D"]
    manifest_path.write_text(json.dumps(manifest))
    assert runner.resume(root)["status"] == "complete"
    assert read(manifest_path)["samples"]["b2AR_and_Nb6B9_vs_Nb80/1"]["observed_chains"] == ["B", "C", "D"]
    assert first_file.stat().st_mtime_ns == first_mtime
    second_file = next(p for p in root.glob("models/*/*/*.pdb") if p != first_file)
    second_mtime = second_file.stat().st_mtime_ns
    first_file.write_text("broken PDB")
    assert runner.resume(root)["status"] == "complete"
    assert first_file.stat().st_size > 100000 and second_file.stat().st_mtime_ns == second_mtime
    assert (root / "ranking.csv").read_bytes() and (root / "heatmap.png").read_bytes().startswith(b"\x89PNG")


def test_chain_identity_in_pdb_and_mmcif(tmp_path):
    cmd = pytest.importorskip("pymol").cmd
    from appraise.nanobody_structure import identify_chains
    source_name = next(iter(panel()["structures"]))
    pdb = tmp_path / source_name
    pdb.write_bytes(gzip.decompress(resource(source_name + ".gz").read_bytes()))
    cmd.reinitialize()
    assert identify_chains(cmd, pdb, "Nb80", "Nb6B9", ["B", "C", "D"]) == ["C", "B", "D"]
    cif = tmp_path / "three_chains.cif"
    cmd.save(str(cif), "inspection")
    cmd.reinitialize()
    assert identify_chains(cmd, cif, "Nb6B9", "Nb80", []) == ["B", "C", "D"]
    cmd.reinitialize()
    with pytest.raises(ValueError, match="Identical"):
        identify_chains(cmd, pdb, "Nb80", "Nb80", [])
    cmd.reinitialize()


def test_unqualified_default_and_missing_gpu_return_actionable_status(tmp_path, monkeypatch):
    result = runner.run(runner.RunConfig(output_dir=str(tmp_path / "auto")))
    assert result["failure_code"] == "unqualified_configuration"
    assert result["expected_predictions"] == 10 and result["completed_predictions"] == 0
    monkeypatch.setattr(runner, "runtime", lambda: {"gpu": None, "packages": {}})
    result = runner.run(runner.RunConfig(output_dir=str(tmp_path / "fresh"), engine="esmfold2-fast"))
    assert result["status"] == "blocked" and result["failure_code"] in {"gpu_unavailable", "pymol_missing"}
    assert "resume" in result["next_action"]
    invalid = subprocess.run([sys.executable, "-m", "appraise", "score", "--predictions-per-input", "0",
                              "--output", str(tmp_path / "invalid")], capture_output=True, text=True)
    assert invalid.returncode == 1 and json.loads(invalid.stdout)["failure_code"] == "configuration_error"


def test_csv_exact_panel_and_freeze_sample_counts(tmp_path):
    path = tmp_path / "panel.csv"
    spec = panel()
    names = ["Nb80", "Nb6B9", "Nb6A9"]
    table = pd.DataFrame([(n, spec["candidates"][n]) for n in names], columns=["peptide_name", "peptide_seq"])
    table.to_csv(path, index=False)
    cfg = runner.RunConfig(example="b2ar-panel", candidates_csv=str(path), engine="esmfold2-fast", seed=17)
    pre = runner.preflight(cfg)
    assert pre["expected_jobs"] == 6 and pre["expected_predictions"] == 30
    samples = [runner._sample_inputs("job", names[:2], pre["settings"], i) for i in range(5)]
    assert [m["seed"] for m in samples] == list(range(17, 22))
    assert {m["sample_index"] for m in samples} == {0}
    table.loc[0, "peptide_seq"] = spec["candidates"]["Nb6B9"]
    table.to_csv(path, index=False)
    with pytest.raises(ValueError, match="mismatch"):
        runner.preflight(cfg)


@pytest.mark.parametrize("engine", ["colabfold", "esmfold2-fast"])
def test_backend_requests_one_identifiable_sample_on_resume(tmp_path, monkeypatch, engine):
    monkeypatch.setattr(nanobody_models.metadata, "version", lambda name: "1.6.3" if name == "colabfold" else "3.1.18")
    monkeypatch.setattr(nanobody_models, "runtime", lambda: {"gpu": {"compute_cap": "7.5"}})
    monkeypatch.setattr(nanobody_models, "_ported_files", lambda path: tmp_path / "entry.py")
    job = "b2AR_and_Nb80_vs_Nb80"
    work = tmp_path / "work"
    def fake_command(cmd, log, env=None):
        if engine == "colabfold":
            assert cmd[cmd.index("--model-order") + 1] == "4"
            assert cmd[cmd.index("--random-seed") + 1] == "17"
            assert cmd[cmd.index("--num-models") + 1] == cmd[cmd.index("--num-seeds") + 1] == "1"
            folder = work / "output"
            model = job + "_unrelaxed_rank_001_alphafold2_multimer_v2_model_4_seed_017.pdb"
            conf = job + "_scores_rank_001_alphafold2_multimer_v2_model_4_seed_017.json"
        else:
            request = read(work / "input.json")
            assert request["modelSeeds"] == [20]
            assert [p["protein"]["id"] for p in request["sequences"]] == [["A", "B"], "C"]
            assert "--use_esm_embeddings" in cmd and "--num_diffusion_samples=1" in cmd
            folder = work / "output" / job.lower() / "seed-20_sample-0"
            model, conf = "native_model.cif", "native_summary_confidences.json"
        folder.mkdir(parents=True)
        (folder / model).write_text("format fixture")
        (folder / conf).write_text(json.dumps({"iptm": 0.5}))
        return {"elapsed_seconds": 1}
    monkeypatch.setattr(nanobody_models, "_checked_command", fake_command)
    model, conf, usage = nanobody_models.predict(engine, job, ["Nb80", "Nb80"], 3, 17, work, tmp_path / "log")
    assert model.is_file() and read(conf)["iptm"] == 0.5 and usage["elapsed_seconds"] == 1
