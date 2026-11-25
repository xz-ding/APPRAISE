"""Candidate modeling routes for the published β₂AR nanobody panel."""

import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
from importlib import metadata
from pathlib import Path
from urllib.request import urlopen

from .nanobody import panel


COLABFOLD_REV = "84c27d9cc500489fd9b97545d2325b9d00f251d5"
PREVIEW_REV = "eb061aa6706959fdffe15e86d27ca08083d6ed8d"
PORT_REV = "4d038d790df99270f81673052aca5aa23b3f9a03"
WEIGHTS_REV = "6527d516857d9006063536809505951e2c30df9f"
SCRIPT_SHA256 = "70a1102908353d75bccfccfa85738f72d053a8fa819554c5e02e6230cc9ba86b"
PROTOCOLS = {
    "saved": {"engine": "AlphaFold2-multimer", "source": "10.22002/kxjgj-tfk18"},
    "colabfold": {"engine": "AlphaFold2-multimer v2", "colabfold": COLABFOLD_REV,
                  "templates": True, "recycles": 3, "msa": "mmseqs2_uniref_env", "relax": False},
    "esmfold2-fast": {"engine": "ColabFold2 preview ESMFold2-Fast", "preview": PREVIEW_REV,
                      "implementation": PORT_REV, "weights_repo": "sokrypton/af3-any-model",
                      "weights_revision": WEIGHTS_REV, "tower": "esmc", "weight_storage": "int8",
                      "attention": "xla", "recycles": 10, "sampling_steps": 15, "samples_per_seed": 1},
}


def runtime():
    result = {"python": sys.version.split()[0], "platform": sys.platform, "gpu": None,
              "colab_release": os.environ.get("COLAB_RELEASE_TAG")}
    if "SC_PHYS_PAGES" in getattr(os, "sysconf_names", {}):
        result["system_memory_bytes"] = os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
    if shutil.which("nvidia-smi"):
        fields = ["name", "memory.free", "memory.total", "compute_cap", "driver_version"]
        probe = subprocess.run(
            ["nvidia-smi", "-i", "0", "--query-gpu=" + ",".join(fields), "--format=csv,noheader,nounits"],
            capture_output=True, text=True, timeout=15,
        )
        if probe.returncode == 0 and probe.stdout.strip():
            result["gpu"] = dict(zip(fields, (field.strip() for field in probe.stdout.splitlines()[0].split(","))))
    result["packages"] = {}
    for name in ("appraise", "pymol-open-source", "colabfold", "alphafold3-colabfold", "jax", "jaxlib",
                 "jax-cuda12-plugin", "jax-cuda13-plugin", "numpy", "pandas", "dm-haiku", "huggingface-hub",
                 "nvidia-cuda-runtime-cu12", "nvidia-cuda-runtime-cu13", "openmm", "colabfold-binaries"):
        try:
            result["packages"][name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            continue
    return result


def _checked_command(command, log, env=None):
    with log.open("a", encoding="utf-8") as stream, tempfile.TemporaryFile(mode="w+") as samples:
        stream.write("Running " + json.dumps(command) + "\n")
        stream.flush()
        monitor = subprocess.Popen(["nvidia-smi", "-i", "0", "--query-gpu=memory.used",
                                    "--format=csv,noheader,nounits", "--loop=1"],
                                   stdout=samples, stderr=subprocess.DEVNULL)
        start = time.monotonic()
        try:
            subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, check=True, env=env)
        except subprocess.CalledProcessError as exc:
            # Include the upstream diagnostic in the structured result as well as the full log.
            with log.open("rb") as errors:
                errors.seek(max(0, errors.seek(0, os.SEEK_END) - 6000))
                lines = errors.read().decode(errors="replace").strip().splitlines()
            if any(marker in " ".join(lines).lower() for marker in ("out of memory", "out of gpu memory", "resource_exhausted")):
                raise MemoryError("The selected GPU ran out of memory during inference. The original settings and completed samples are preserved.") from exc
            raise RuntimeError(lines[-1] if lines else str(exc)) from exc
        finally:
            monitor.terminate()
            monitor.wait(timeout=10)
            samples.seek(0)
            values = [int(line.strip()) for line in samples if line.strip().isdigit()]
            usage = {"elapsed_seconds": round(time.monotonic() - start, 3),
                     "peak_observed_device_mib": max(values, default=None), "sampling_interval_seconds": 1}
            stream.write(json.dumps(usage) + "\n")
    return usage


def _ported_files(cache):
    from huggingface_hub import hf_hub_download
    from alphafold3.constants import ccd_fetch

    package = Path(metadata.distribution("alphafold3-colabfold").locate_file("alphafold3"))
    converters = package / "constants" / "converters"
    ccd, sets = converters / "ccd.pickle", converters / "chemical_component_sets.pickle"
    if not (ccd.is_file() and sets.is_file()):
        converters.mkdir(parents=True, exist_ok=True)
        ccd_fetch.write_pickles(ccd_fetch.codes_for_input(extra=[]), str(ccd), str(sets),
                                libcifpp_dir=str(package.parent / "share" / "libcifpp"))

    cache.mkdir(parents=True, exist_ok=True)
    entry = cache / "run_alphafold.py"
    if not entry.is_file() or hashlib.sha256(entry.read_bytes()).hexdigest() != SCRIPT_SHA256:
        url = f"https://raw.githubusercontent.com/sokrypton/alphafold3/{PORT_REV}/run_alphafold.py"
        with urlopen(url, timeout=60) as response:
            source = response.read()
        if hashlib.sha256(source).hexdigest() != SCRIPT_SHA256:
            raise ValueError("The upstream ColabFold2 entry point failed its pinned checksum.")
        entry.write_bytes(source)
    assets = {
        "esmfold2/esmfold2_fast.int8.bin.zst": "esmfold2_fast-int8",
        "esmfold2/esmfold2_fast.lm.npz": "esmfold2_fast-int8",
        "lm/esmc.bin.zst": "esmc",
    }
    for remote, folder in assets.items():
        local = Path(hf_hub_download("sokrypton/af3-any-model", remote, revision=WEIGHTS_REV))
        dest = cache / folder / Path(remote).name
        dest.parent.mkdir(exist_ok=True)
        if not dest.is_symlink() or dest.resolve() != local.resolve():
            if dest.exists() or dest.is_symlink():
                raise ValueError(f"Model cache has an unexpected existing file at {dest}.")
            dest.symlink_to(local.resolve())
    return entry


def predict(engine, job, names, index, base_seed, work, log):
    """Produce exactly one identifiable structure for a fixed published competition."""
    spec = panel()
    seqs = [spec["candidates"][name] for name in names] + [spec["receptor"]["sequence"]]
    path = Path(work)
    path.mkdir(parents=True, exist_ok=True)
    if engine == "colabfold":
        if metadata.version("colabfold") != "1.6.3":
            raise ValueError("The candidate AlphaFold2 route requires ColabFold 1.6.3.")
        fasta = path / f"{job}.fasta"
        fasta.write_text(f">{job}\n" + ":".join(seqs) + "\n")
        out = path / "output"
        command = [sys.executable, "-m", "colabfold.batch", str(fasta), str(out),
                   "--model-type", "alphafold2_multimer_v2", "--templates", "--num-models", "1",
                   "--model-order", str(index + 1), "--num-seeds", "1", "--random-seed", str(base_seed),
                   "--num-recycle", "3", "--msa-mode", "mmseqs2_uniref_env"]
        usage = _checked_command(command, log)
        products = list(out.glob(f"*_unrelaxed_rank_*_model_{index+1}_seed_*.pdb"))
        cofactors = list(out.glob("*_scores_rank_*.json"))
    elif engine == "esmfold2-fast":
        if metadata.version("alphafold3-colabfold") != "3.1.18":
            raise ValueError("The candidate ColabFold2 route requires alphafold3-colabfold 3.1.18.")
        cache = Path.home() / ".cache" / "appraise" / WEIGHTS_REV
        started = time.monotonic()
        entry = _ported_files(cache)
        setup_seconds = time.monotonic() - started
        seed = base_seed + index
        # AF3 JSON expresses identical copies as an array of assigned chain IDs.
        groups = {}
        for chain_id, seq in zip("ABC", seqs):
            groups.setdefault(seq, []).append(chain_id)
        proteins = [{"protein": {"id": ids[0] if len(ids) == 1 else ids, "sequence": seq,
                                 "unpairedMsa": f">query\n{seq}\n", "pairedMsa": "", "templates": []}}
                    for seq, ids in groups.items()]
        query = {"name": job.lower(), "dialect": "alphafold3", "version": 1, "modelSeeds": [seed], "sequences": proteins}
        inp = path / "input.json"
        inp.write_text(json.dumps(query))
        out = path / "output"
        command = [sys.executable, str(entry), f"--json_path={inp}", f"--output_dir={out}",
                   "--model=esmfold2_fast", "--use_esm_embeddings", "--norun_data_pipeline",
                   f"--model_dir={cache / 'esmfold2_fast-int8'}", "--nodownload_weights", "--weights_precision=int8",
                   "--num_recycles=10", "--num_sampling_steps=15", "--num_diffusion_samples=1",
                   "--flash_attention_implementation=xla", "--force_output_dir", "--nocompress_large_output_files",
                   f"--cache_dir={cache / 'compiled'}"]
        env = dict(os.environ, AF3_WEIGHTS_DIR=str(cache))
        cap = float(runtime()["gpu"]["compute_cap"])
        flag = "--xla_disable_hlo_passes=custom-kernel-fusion-rewriter" if cap < 8 else "--xla_gpu_enable_triton_gemm=false"
        env["XLA_FLAGS"] = (env.get("XLA_FLAGS", "") + " " + flag).strip()
        usage = _checked_command(command, log, env)
        usage["weight_and_script_setup_seconds"] = round(setup_seconds, 3)
        sample = out / job.lower() / f"seed-{seed}_sample-0"
        products = list(sample.glob("*.cif"))
        cofactors = list(sample.glob("*_summary_confidences.json"))
    else:
        raise ValueError(f"Unknown nanobody modeling engine {engine!r}.")
    if len(products) != 1 or len(cofactors) != 1:
        raise ValueError(f"Expected one structure and confidence record for {job}, prediction {index + 1}.")
    confidence = json.loads(cofactors[0].read_text())
    iptm = confidence.get("iptm")
    if not isinstance(iptm, (float, int)) or not math.isfinite(iptm) or not 0 <= iptm <= 1:
        raise ValueError(f"Missing native multichain confidence in {cofactors[0]}.")
    return products[0], cofactors[0], usage
