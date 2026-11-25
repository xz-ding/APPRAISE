"""Run and resume the published β₂AR nanobody competitions."""

import contextlib
import csv
import fcntl
import gzip
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import time
import uuid
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from . import __version__
from .nanobody import comparison, evaluate, ordered_inputs, panel, reference_measurements, resource
from .nanobody_models import PROTOCOLS, predict, runtime


@dataclass(frozen=True)
class RunConfig:
    example: str = "b2ar-two"
    engine: str = "auto"
    output_dir: str = "runs/b2ar-demo"
    predictions_per_input: int = 5
    seed: int = 0
    candidates_csv: str | None = None


class RunIssue(Exception):
    def __init__(self, code, message, action, status="failed"):
        super().__init__(message)
        self.code, self.action, self.status = code, action, status


def _hash(content):
    return hashlib.sha256(content).hexdigest()


def _revision(*modules):
    return _hash(b"".join(Path(__file__).with_name(name + ".py").read_bytes() for name in modules))


def _json(path, data):
    _write(path, (json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n").encode())


def _write(path, content):
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
        tmp = Path(stream.name)
        try:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise
    os.replace(tmp, path)


def _read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def _current(path, digest):
    return bool(digest) and path.is_file() and _hash(path.read_bytes()) == digest


def _now():
    return datetime.now(timezone.utc).isoformat()


def _configuration(config):
    spec = panel()
    known = spec["candidates"]
    if config.example not in {"b2ar-two", "b2ar-saved", "b2ar-reference", "b2ar-six", "b2ar-panel"}:
        raise ValueError("Choose b2ar-two, b2ar-saved, b2ar-reference, b2ar-six or b2ar-panel.")
    if type(config.predictions_per_input) is not int or not 1 <= config.predictions_per_input <= 5:
        raise ValueError("predictions_per_input must be an integer from 1 through 5.")
    if type(config.seed) is not int or not 0 <= config.seed <= 2**31 - 5:
        raise ValueError("seed must be a nonnegative 31-bit integer with room for five samples.")
    archived = config.example == "b2ar-reference"
    full = config.example == "b2ar-six" or archived
    if full and config.predictions_per_input != 5:
        raise ValueError("The six-nanobody reference and benchmark require five predictions per input.")
    names = list(known) if full else ["Nb6B9", "Nb80"]
    if config.example == "b2ar-panel":
        if not config.candidates_csv:
            raise ValueError("Supply --candidates-csv with peptide_name and peptide_seq columns.")
        with Path(config.candidates_csv).open(encoding="utf-8-sig", newline="") as stream:
            rows = list(csv.DictReader(stream))
        names = [row.get("peptide_name", "").strip() for row in rows]
        if len(names) < 2 or len(set(names)) != len(names):
            raise ValueError("Select at least two distinct published β₂AR nanobodies.")
        for name, row in zip(names, rows):
            if name not in known or row.get("peptide_seq", "").strip().upper() != known[name]:
                raise ValueError(f"Name or exact published β₂AR sequence mismatch for {name!r}.")
    elif config.candidates_csv:
        raise ValueError("Choose b2ar-panel to use a candidates CSV.")
    if config.engine not in {"auto", "saved", "reference", "colabfold", "esmfold2-fast"}:
        raise ValueError("Choose auto, colabfold, esmfold2-fast, saved or reference.")
    engine = config.engine
    if archived or config.example == "b2ar-saved":
        required = "reference" if archived else "saved"
        if engine not in {"auto", required} or config.seed != 0:
            raise ValueError("Archived inputs have a fixed producer and sample identities. Keep the preset engine and seed 0.")
        engine = required
    elif engine in {"saved", "reference"}:
        raise ValueError("Use b2ar-saved for archived structures or b2ar-reference for archived measurements.")
    protocol = PROTOCOLS.get(engine, {"engine": engine, "source": spec["source"]})
    settings = {"example": config.example, "engine": engine, "names": names, "self_pairs": full,
                "sequences": {name: known[name] for name in names}, "receptor": spec["receptor"],
                "predictions_per_input": config.predictions_per_input, "seed": config.seed, "protocol": protocol}
    return settings


def _plan(settings):
    jobs = ordered_inputs(settings["names"], settings["self_pairs"])
    count = len(jobs) * settings["predictions_per_input"]
    return jobs, {"expected_jobs": len(jobs), "expected_predictions": count, "expected_measurements": count * 2}


def preflight(config):
    settings = _configuration(config)
    _, counts = _plan(settings)
    hardware = runtime()
    return dict(counts, settings=settings, runtime=hardware, output_dir=str(Path(config.output_dir).resolve()),
                selection_reason=_reason(settings["engine"]))


def _reason(engine):
    return {
        "auto": "No fresh-folding configuration has completed free-Colab qualification. Choose a validation engine explicitly or score the saved structures.",
        "reference": "Re-score the published six-nanobody measurements on CPU.",
        "saved": "Re-measure the ten published reciprocal structures on CPU.",
        "colabfold": "Explicit AlphaFold2-multimer v2 validation run with templates.",
        "esmfold2-fast": "Explicit experimental ColabFold2 preview ESMFold2-Fast validation run.",
    }[engine]


def _sample_inputs(job, names, settings, index):
    engine = settings["engine"]
    receptor = settings["receptor"]
    chains = [{"role": role, "input_id": chain, "name": name, "sequence_sha256": _hash(seq.encode())}
              for role, chain, name, seq in zip(("candidate_1", "candidate_2", "receptor"), "ABC",
                  [*names, receptor["name"]], [*(settings["sequences"][n] for n in names), receptor["sequence"]])]
    seed = None if engine == "saved" else settings["seed"]
    if engine == "esmfold2-fast":
        seed += index
    return {"job": job, "chains": chains, "checkpoint_model": index + 1 if engine == "colabfold" else None,
            "seed": seed,
            "sample_index": 0 if engine == "esmfold2-fast" else index,
            "ordinal": index + 1, "engine": engine, "protocol": settings["protocol"]}


def _worker(args, log):
    with log.open("a", encoding="utf-8") as stream:
        completed = subprocess.run([sys.executable, "-m", "appraise.nanobody_structure", *map(str, args)],
                                   stdout=stream, stderr=subprocess.STDOUT)
    if completed.returncode:
        tail = log.read_bytes()[-3500:].decode(errors="replace").strip().splitlines()
        raise ValueError(tail[-1] if tail else f"PyMOL exited with code {completed.returncode}.")


def _prepare_structure(folder, scratch, job, names, index, settings, record, log):
    engine = settings["engine"]
    target = folder / (f"{job}_unrelaxed_model_{index + 1}" + (".cif" if engine == "esmfold2-fast" else ".pdb"))
    confidence = folder / "confidence.json"
    if (_current(target, record.get("structure_sha256")) and record.get("observed_chains")
            and (engine == "saved" or _current(confidence, record.get("confidence_sha256")))):
        return target
    if engine == "saved":
        structures = panel()["structures"]
        archived = sorted(n for n in structures if n.startswith(job + "_unrelaxed_"))[index]
        data = gzip.decompress(resource(archived + ".gz").read_bytes())
        if _hash(data) != structures[archived]:
            raise ValueError(f"The bundled structure failed its archive checksum: {archived}")
        record["source_file"] = archived
        declared = "B,C,D"
    else:
        if not runtime()["gpu"]:
            raise RunIssue("gpu_unavailable", "New predictions require an attached NVIDIA GPU.",
                           "Connect a free Colab GPU and resume, or start b2ar-saved in a separate directory.", "blocked")
        distribution = "colabfold" if engine == "colabfold" else "alphafold3-colabfold"
        if distribution not in runtime()["packages"]:
            raise RunIssue("engine_missing", f"The selected engine requires {distribution}.",
                           "Run the candidate notebook setup for the selected engine, then resume.", "blocked")
        model, scores, usage = predict(engine, job, names, index, settings["seed"], scratch, log)
        record["resources"] = usage
        data = model.read_bytes()
        _write(confidence, scores.read_bytes())
        record["confidence_sha256"] = _hash(confidence.read_bytes())
        declared = "A,B,C"
    # Inspect before publishing an inference checkpoint. Reuse is tied to these exact bytes.
    candidate = scratch / target.name
    candidate.write_bytes(data)
    identity = scratch / "chains.json"
    _worker(["--inspect", candidate, identity, *names, declared], log)
    _write(target, data)
    record.update(structure_sha256=_hash(data), observed_chains=_read(identity), structure=str(target.relative_to(folder.parents[2])),
                  plddt_source="native model confidence, 0-100", structure_verified_at=_now())
    record.pop("measurement_sha256", None)
    return target


def _report(root, frame, settings):
    # The scorer prints historical diagnostics. Keep CLI stdout exclusively machine-readable.
    with (root / "run.log").open("a", encoding="utf-8") as log, contextlib.redirect_stdout(log):
        scored, matrix, ranking = evaluate(frame, settings["names"], settings["predictions_per_input"],
                                          self_pairs=settings["self_pairs"], historical=settings["engine"] == "reference")
    order = ranking["peptide_name"].tolist()
    grid = matrix.reindex(index=order, columns=order)
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from io import BytesIO
    import numpy as np
    fig = Figure(figsize=(max(4.8, 0.72 * len(order) + 2), max(4, 0.68 * len(order) + 1.5)), layout="constrained")
    FigureCanvasAgg(fig)
    ax = fig.subplots()
    limit = max(1.0, np.nanmax(np.abs(grid.to_numpy())))
    from matplotlib import colormaps
    cmap = colormaps["RdBu_r"].with_extremes(bad="#e4e4e4")
    im = ax.imshow(grid.to_numpy(), cmap=cmap, vmin=-limit, vmax=limit)
    ax.set(xticks=range(len(order)), xticklabels=order, yticks=range(len(order)), yticklabels=order,
           xlabel="Competitor", ylabel="Candidate", title="β₂AR nanobody competition")
    ax.tick_params(axis="x", labelrotation=35)
    for i, j in zip(*np.where(grid.isna())):
        ax.text(j, i, "—", ha="center", va="center")
    fig.colorbar(im, ax=ax, label="Mean APPRAISE ΔB")
    buffer = BytesIO()
    fig.savefig(buffer, format="png", dpi=140)
    _write(root / "heatmap.png", buffer.getvalue())
    for data, filename in [(scored, "measurements.csv"), (ranking, "ranking.csv")]:
        _write(root / filename, data.to_csv(index=False).encode())
    return comparison(ranking) if settings["example"] in {"b2ar-reference", "b2ar-six"} else None


def _execute(root, settings, manifest, result, report_progress):
    jobs, totals = _plan(settings)
    engine = settings["engine"]
    manifest_path = root / "run_manifest.json"
    log = root / "run.log"
    def checkpoint(stage, active_input=None):
        result.update(stage=stage, active_input=active_input)
        with log.open("a", encoding="utf-8") as stream:
            stream.write(f"{_now()} {stage} {active_input or ''} {result['completed_predictions']}/{totals['expected_predictions']}\n")
        _json(manifest_path, manifest)
        _json(root / "result.json", result)
        if report_progress:
            report_progress(result.copy())

    checkpoint("preflight")
    if engine == "auto":
        raise RunIssue("unqualified_configuration", _reason(engine),
                       "Use b2ar-saved on CPU, or select a validation engine in a new run directory.", "blocked")
    if engine == "reference":
        frame = reference_measurements()
        result.update(completed_jobs=len(jobs), completed_predictions=frame["model_name"].nunique(), completed_measurements=len(frame))
    else:
        if importlib.util.find_spec("pymol") is None:
            raise RunIssue("pymol_missing", "Headless PyMOL is required to verify and measure structures.",
                           "Install this APPRAISE revision with the structures extra, then resume.", "blocked")
        chunks = []
        quant_revision = manifest["quantification_revision"]
        for job, names in jobs.items():
            for i in range(settings["predictions_per_input"]):
                key = f"{job}/{i + 1}"
                identity = _sample_inputs(job, names, settings, i)
                fingerprint = _hash(json.dumps([identity, manifest["inference_revision"]], sort_keys=True).encode())
                sample = manifest["samples"].setdefault(key, {})
                if sample.get("fingerprint") != fingerprint:
                    sample.clear()
                    sample.update(identity=identity, fingerprint=fingerprint)
                folder = root / "models" / job / str(i + 1)
                folder.mkdir(parents=True, exist_ok=True)
                checkpoint("structure", key)
                with tempfile.TemporaryDirectory(prefix="appraise-model-") as scratch:
                    started = time.monotonic()
                    saved_at = sample.get("structure_verified_at")
                    structure = _prepare_structure(folder, Path(scratch), job, names, i, settings, sample, log)
                    if sample.get("structure_verified_at") != saved_at:
                        sample["structure_seconds"] = time.monotonic() - started
                    checkpoint("quantification", key)
                    measurements = folder / "measurements.csv"
                    if sample.get("quantification_revision") != quant_revision or not _current(measurements, sample.get("measurement_sha256")):
                        with tempfile.TemporaryDirectory(prefix="appraise-quant-") as quant:
                            rows = Path(quant) / "rows.csv"
                            _worker([structure, rows, *names, ",".join(sample["observed_chains"])], log)
                            _write(measurements, rows.read_bytes())
                        sample.update(quantification_revision=quant_revision, measurement_sha256=_hash(measurements.read_bytes()))
                    chunks.append(pd.read_csv(measurements))
                result["completed_predictions"] += 1
                result["completed_measurements"] += len(chunks[-1])
                sample["complete"] = True
                if i + 1 == settings["predictions_per_input"]:
                    result["completed_jobs"] += 1
                checkpoint("measured", key)
        frame = pd.concat(chunks, ignore_index=True)
    checkpoint("ranking")
    scientific = _report(root, frame, settings)
    if scientific is not None:
        scientific["evidence"] = "published AlphaFold2 measurements" if engine == "reference" else "newly generated structures"
        result["comparison"] = scientific
    result["artifacts"] = {name: {"path": str(root / name), "sha256": _hash((root / name).read_bytes())}
                           for name in ("measurements.csv", "ranking.csv", "heatmap.png")}
    checkpoint("finished")


def run(config, *, progress=None):
    """Execute a run, or verify and resume it when the same directory already exists."""
    root = Path(config.output_dir).expanduser().resolve()
    result = {"schema_version": 1, "run_id": None, "status": "failed", "selected_engine": config.engine,
              "stage": "configuration", "completed_jobs": 0, "completed_predictions": 0, "completed_measurements": 0,
              "paths": {name: str(root / name) for name in ("run_manifest.json", "result.json", "run.log")}}
    manifest = None
    started = time.monotonic()
    lock = None
    try:
        settings = _configuration(config)
        _, counts = _plan(settings)
        result.update(counts, selected_engine=settings["engine"], selection_reason=_reason(settings["engine"]))
        root.mkdir(parents=True, exist_ok=True)
        lock = (root / ".appraise.lock").open("a")
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RunIssue("run_busy", "Another process is writing this run directory.", "Wait for that process, then resume.", "blocked")
        manifest_path = root / "run_manifest.json"
        if manifest_path.exists():
            existing = _read(manifest_path)
            if (existing.get("schema_version") != 1 or not isinstance(existing.get("samples"), dict)
                    or not isinstance(existing.get("segments"), list) or not existing.get("run_id")):
                raise RunIssue("invalid_manifest", "The saved run manifest is unreadable or has an unsupported schema.",
                               "Restore its last durable copy, or use a new directory.")
            if existing["settings"] != settings:
                raise RunIssue("config_mismatch", "The saved run uses different inputs or settings.",
                               "Resume the existing run or choose a new output directory.")
            manifest = existing
        else:
            if any(p.name != ".appraise.lock" for p in root.iterdir()):
                raise ValueError("The output directory already contains other files. Choose a new directory.")
            manifest = {"schema_version": 1, "run_id": str(uuid.uuid4()), "settings": settings,
                        "config": dict(asdict(config), output_dir=str(root)), "samples": {}, "segments": []}
            if config.candidates_csv:
                normalized = pd.DataFrame(list(settings["sequences"].items()), columns=["peptide_name", "peptide_seq"])
                _write(root / "candidates.csv", normalized.to_csv(index=False).encode())
                manifest["config"]["candidates_csv"] = "candidates.csv"
        hardware = runtime()
        # Keep validation/quantification edits independent of expensive inference.
        manifest.update(inference_revision="b2ar-structure-v1/" + _revision("nanobody_models"),
                        quantification_revision=_revision("pymol_quantify_peptide_binding", "nanobody_structure") + "/" + hardware["packages"].get("pymol-open-source", "unknown"),
                        scoring_revision=_revision("nanobody", "score_calculation", "utilities"))
        segment = {"start": _now(), "environment": hardware, "package_version": __version__}
        manifest["segments"].append(segment)
        result.update(run_id=manifest["run_id"], status="interrupted")
        _execute(root, settings, manifest, result, progress)
        result.update(status="complete", failure_code=None, next_action=None)
    except KeyboardInterrupt:
        result.update(status="interrupted", failure_code="interrupted", error="Run interrupted.", next_action="Resume the same run directory.")
    except RunIssue as exc:
        result.update(status=exc.status, failure_code=exc.code, error=str(exc), next_action=exc.action)
    except (OSError, ValueError, KeyError, ImportError, subprocess.SubprocessError, RuntimeError) as exc:
        stage = result["stage"]
        action = "Check the inputs and output directory." if stage == "configuration" else "Read run.log, correct the reported error and resume the same directory."
        result.update(status="failed", failure_code=f"{stage}_error", error=str(exc), next_action=action)
    finally:
        try:
            result["elapsed_seconds"] = round(time.monotonic() - started, 3)
            if manifest is not None:
                if manifest["segments"]:
                    manifest["segments"][-1].update(end=_now(), status=result["status"])
                manifest["last_status"] = result["status"]
                _json(root / "result.json", result)
                _json(root / "run_manifest.json", manifest)
        except OSError as exc:
            result.update(status="failed", failure_code="persistence_error", error=str(exc),
                          next_action="Restore write access to the destination and resume.")
        finally:
            if lock is not None:
                lock.close()
    return result


def resume(run_dir, *, progress=None):
    root = Path(run_dir).expanduser().resolve()
    try:
        config = _read(root / "run_manifest.json")["config"]
        config["output_dir"] = str(root)
        if config.get("candidates_csv"):
            config["candidates_csv"] = str(root / "candidates.csv")
        return run(RunConfig(**config), progress=progress)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        return {"schema_version": 1, "status": "failed", "failure_code": "invalid_manifest",
                "error": str(exc), "next_action": "Choose the directory containing the original run_manifest.json."}
