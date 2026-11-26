"""Verify the built wheel and sdist against the package source, then exercise the installed wheel."""

import argparse
import hashlib
import json
import sys
import tarfile
import zipfile
from pathlib import Path
from tempfile import TemporaryDirectory


def verify_files(source, dist):
    package = source / "appraise"
    # The legacy misc_utilities directory holds standalone source-checkout tools.
    shipped = [*package.glob("*.py"), *(package / "examples/b2ar").glob("*")]
    expected = {p.relative_to(source).as_posix(): p.read_bytes()
                for p in shipped if p.is_file() and not p.name.startswith(".")}
    wheel, = dist.glob("*.whl")
    sdist, = dist.glob("*.tar.gz")
    with zipfile.ZipFile(wheel) as archive:
        actual = {name: archive.read(name) for name in archive.namelist() if name.startswith("appraise/") and not name.endswith("/")}
        assert actual == expected, "Wheel package files differ from the candidate source."
    with tarfile.open(sdist) as archive:
        actual = {item.name.split("/", 1)[1]: archive.extractfile(item).read() for item in archive
                  if item.isfile() and "/appraise/" in item.name}
        assert actual == expected, "Source distribution package files differ from the candidate source."
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (wheel, sdist)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--dist", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    import appraise
    from appraise.runner import RunConfig, resume, run
    from importlib.metadata import version
    assert Path(appraise.__file__).resolve().is_relative_to(Path(sys.prefix).resolve())
    assert appraise.__version__ == version("appraise")
    artifacts = verify_files(args.source.resolve(), args.dist.resolve())
    report = {"version": appraise.__version__, "artifacts": artifacts, "examples": {}}
    with TemporaryDirectory() as tmp:
        for example, expected in (("b2ar-reference", (36, 180, 360)), ("b2ar-saved", (2, 10, 20))):
            root = Path(tmp) / example
            result = run(RunConfig(example=example, output_dir=str(root)))
            assert result["status"] == "complete", result
            assert tuple(result[key] for key in ("completed_jobs", "completed_predictions", "completed_measurements")) == expected
            from csv import DictReader
            with (root / "ranking.csv").open() as f:
                order = [r["peptide_name"] for r in DictReader(f)]
            assert order == (["Nb6B9", "Nb6A5", "Nb6A9", "Nb6A10", "Nb6B11", "Nb80"] if example == "b2ar-reference" else ["Nb6B9", "Nb80"])
            if example == "b2ar-reference":
                assert abs(result["comparison"]["spearman_rho"] - 31 / 35) < 1e-12
            completed = resume(root)
            assert completed["status"] == "complete" and completed["artifacts"] == result["artifacts"]
            report["examples"][example] = {
                "jobs": expected[0], "structures": expected[1], "measurements": expected[2], "order": order,
                "comparison": result.get("comparison"), "ranking_sha256": result["artifacts"]["ranking.csv"]["sha256"],
            }
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
