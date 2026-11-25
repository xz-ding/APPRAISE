"""Command-line interface to the shared APPRAISE nanobody workflow."""

import argparse
import json
import sys

from .runner import RunConfig, resume, run


def main(argv=None):
    parser = argparse.ArgumentParser(prog="python -m appraise")
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("run", "score"):
        cmd = commands.add_parser(name, help="Model published competitions" if name == "run" else "Analyze archived evidence on CPU")
        cmd.add_argument("--example", default="b2ar-two" if name == "run" else "b2ar-saved",
                         choices=["b2ar-two", "b2ar-six", "b2ar-panel"] if name == "run" else ["b2ar-saved", "b2ar-reference"])
        cmd.add_argument("--output", required=True, dest="output_dir")
        cmd.add_argument("--predictions-per-input", default=5, type=int)
        if name == "run":
            cmd.add_argument("--engine", default="auto", choices=["auto", "colabfold", "esmfold2-fast"])
            cmd.add_argument("--seed", default=0, type=int)
            cmd.add_argument("--candidates-csv")
    commands.add_parser("resume", help="Verify and resume a run").add_argument("--run-dir", required=True)
    args = vars(parser.parse_args(argv))
    command = args.pop("command")
    def progress(event):
        print(f"{event['stage']} {event.get('active_input') or ''} "
              f"{event['completed_predictions']}/{event['expected_predictions']} · {event['elapsed_seconds']:.0f}s",
              file=sys.stderr, flush=True)
    result = resume(args["run_dir"], progress=progress) if command == "resume" else run(RunConfig(**args), progress=progress)
    print(json.dumps(result, indent=2, allow_nan=False))
    return {"complete": 0, "failed": 1, "blocked": 3, "interrupted": 130}[result["status"]]


if __name__ == "__main__":
    sys.exit(main())
