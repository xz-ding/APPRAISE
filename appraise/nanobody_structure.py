"""Quantify one published β₂AR competition in an isolated, headless PyMOL process."""

import argparse
import csv
import json
import shutil
import tempfile
from pathlib import Path

import numpy as np

from .nanobody import panel


def identify_chains(cmd, path, left, right, declared):
    spec = panel()
    expected = [spec["candidates"][left], spec["candidates"][right], spec["receptor"]["sequence"]]
    cmd.load(str(path), "inspection")
    if cmd.count_states("inspection") != 1:
        raise ValueError("Expected one model per saved structure.")
    observed = {}
    for chain in cmd.get_chains("inspection"):
        if not chain.isalnum():
            raise ValueError(f"Unsupported chain identifier {chain!r}.")
        selection = f"inspection and chain {chain} and polymer.protein"
        observed[chain] = "".join(cmd.get_fastastr(selection).splitlines()[1:])
        coords = cmd.get_coords(selection)
        if coords is None or not np.isfinite(coords).all():
            raise ValueError(f"Missing or nonfinite coordinates for chain {chain}.")
        atoms = cmd.get_model(selection + " and name CA").atom
        if len(atoms) != len(observed[chain]):
            raise ValueError(f"Residue coordinates are incomplete in chain {chain}.")
        confidence = np.array([atom.b for atom in atoms])
        if not np.isfinite(confidence).all() or (confidence < 0).any() or (confidence > 100).any():
            raise ValueError(f"Invalid native pLDDT on the 0–100 scale for chain {chain}.")
    if len(observed) != 3:
        raise ValueError(f"Expected three protein chains in {path.name}, found {list(observed)}.")
    if declared and [observed.get(k) for k in declared] == expected:
        return list(declared)
    if len(set(expected)) < len(expected):
        raise ValueError("Identical nanobody copies require matching chain IDs from the producer.")
    by_sequence = {sequence: chain for chain, sequence in observed.items()}
    if set(by_sequence) != set(expected):
        raise ValueError("Written chain sequences differ from the manuscript inputs.")
    return [by_sequence[s] for s in expected]


def measure(path, destination, left, right, declared):
    from pymol import cmd
    from .pymol_quantify_peptide_binding import quantify_results_folder

    chains = identify_chains(cmd, path, left, right, declared)
    cmd.delete("all")
    # Give the folder-based quantifier exactly one input, even if a checkpoint holds other files.
    with tempfile.TemporaryDirectory() as model_dir:
        shutil.copyfile(path, Path(model_dir) / path.name)
        quantify_results_folder(
            model_dir, output_path=str(destination.parent), time_stamp=False,
            receptor_chain=chains[-1], chain_order=chains, anchor_site="N-term", glycine_linkers=False,
        )
    produced = destination.parent / "database_APPRAISE_measurements.csv"
    with produced.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 2 or [r["peptide_seq"] for r in rows] != [panel()["candidates"][n] for n in (left, right)]:
        raise ValueError("PyMOL did not produce both expected full-length candidate perspectives.")
    for row in rows:
        row["model_name"] = path.stem
    with destination.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    return chains


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inspect", action="store_true")
    parser.add_argument("structure", type=Path)
    parser.add_argument("measurements", type=Path)
    parser.add_argument("left")
    parser.add_argument("right")
    parser.add_argument("declared_chains")
    args = parser.parse_args()
    declared = args.declared_chains.split(",")
    if args.inspect:
        from pymol import cmd
        chains = identify_chains(cmd, args.structure, args.left, args.right, declared)
        args.measurements.write_text(json.dumps(chains))
    else:
        measure(args.structure, args.measurements, args.left, args.right, declared)
