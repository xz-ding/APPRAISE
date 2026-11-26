import json
import subprocess
import sys
import time
from pathlib import Path
from unittest.mock import patch


BOOK = Path(__file__).resolve().parents[1] / "notebooks/Colab_APPRAISE_candidate.ipynb"


def source(cell_id, book=BOOK):
    cell, = (c for c in json.loads(book.read_text())["cells"] if c["id"] == cell_id)
    return "".join(cell["source"])


def test_cpu_setup_skips_all_folding_installs(tmp_path):
    namespace = dict(engine="reference", fresh=False, project=tmp_path, workflow="Published six-nanobody measurements",
                     Path=Path, subprocess=subprocess, sys=sys, time=time, json=json)
    def receipt(*args, **kwargs):
        return json.dumps({"source": namespace["revision"], "version": "1.3.2", "python": "3.12"})
    with patch("subprocess.run", return_value=subprocess.CompletedProcess([], 0)) as calls, patch("subprocess.check_output", receipt):
        exec(compile(source("install"), "install", "exec"), namespace)
    args = [str(arg) for call in calls.call_args_list for arg in call.args[0]]
    assert not any("alphafold3-colabfold" in arg or "sokrypton" in arg or "pymol" in arg for arg in args)
    assert any(f"@{namespace['revision']}" in arg for arg in args)
    assert len(namespace["revision"]) == 40
    assert json.loads(next(tmp_path.rglob("installation.json")).read_text())["source"] == namespace["revision"]


def test_run_cell_invokes_real_shared_command_and_resumes(tmp_path):
    namespace = dict(python=sys.executable, run_dir=tmp_path / "notebook", example="b2ar-reference",
                     fresh=False, subprocess=subprocess, json=json)
    code = compile(source("run"), "run", "exec")
    exec(code, namespace)
    first = namespace["result"]
    assert first["status"] == "complete" and first["completed_measurements"] == 360
    exec(code, namespace)
    assert "resume" in namespace["command"]
    assert namespace["result"]["run_id"] == first["run_id"]
    assert namespace["result"]["artifacts"] == first["artifacts"]


def test_public_notebook_installs_the_exact_pypi_release_and_runs_saved_example(tmp_path):
    book = BOOK.with_name("Colab_APPRAISE.ipynb")
    project = tmp_path / "public"
    namespace = dict(engine="saved", fresh=False, project=project, workflow="Saved structures on CPU",
                     Path=Path, subprocess=subprocess, sys=sys, time=time, json=json)
    receipt = json.dumps({"version": "1.3.2", "source": "PyPI", "python": "3.12"})
    with patch("subprocess.run", return_value=subprocess.CompletedProcess([], 0)) as calls, patch(
        "subprocess.check_output", return_value=receipt
    ):
        exec(compile(source("install", book), "public install", "exec"), namespace)
    requirements = [cmd for cmd in calls.call_args_list if "appraise[structures]==1.3.2" in cmd.args[0]]
    assert len(requirements) == 1
    assert "https://pypi.org/simple" in requirements[0].args[0]
    assert json.loads(next(project.rglob("installation.json")).read_text())["source"] == "PyPI"

    namespace.update(python=sys.executable, run_dir=tmp_path / "saved", example="b2ar-saved")
    exec(compile(source("run", book), "public run", "exec"), namespace)
    result = namespace["result"]
    assert result["status"] == "complete" and result["completed_measurements"] == 20
    assert {"ranking.csv", "heatmap.png"} <= set(result["artifacts"])
