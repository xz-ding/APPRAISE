"""
Utilities for preparing APPRAISE-compatible modeling inputs and structure files.

These helpers keep the APPRAISE job name as the stable identifier across
different protein-structure backends so that modeled structures can be routed
back into the existing APPRAISE quantification and scoring pipeline.
"""

from __future__ import annotations

import glob
import json
import os
import re
from pathlib import Path
from string import ascii_uppercase, ascii_lowercase


STRUCTURE_FILE_EXTENSIONS = (".pdb", ".cif", ".mmcif")
MULTICHAIN_BACKENDS = {"boltz1", "chai1", "openfold3"}
_CHAIN_IDS = list(ascii_uppercase + ascii_lowercase)
_AUXILIARY_DIRECTORY_PREFIXES = ("templates_",)
_AUXILIARY_DIRECTORY_SUFFIXES = ("_env", "_pairgreedy")
_BACKEND_ALIASES = {
    "af2multimer": "legacy_fasta",
    "af2multimerv1": "legacy_fasta",
    "af2multimerv2": "legacy_fasta",
    "af2multimerv3": "legacy_fasta",
    "alphafold2multimer": "legacy_fasta",
    "alphafold2multimerv1": "legacy_fasta",
    "alphafold2multimerv2": "legacy_fasta",
    "alphafold2multimerv3": "legacy_fasta",
    "alphafoldmultimer": "legacy_fasta",
    "colabfold": "legacy_fasta",
    "esmfold": "legacy_fasta",
    "esmfold2": "legacy_fasta",
    "legacy": "legacy_fasta",
    "legacyfasta": "legacy_fasta",
    "boltz": "boltz1",
    "boltz1": "boltz1",
    "chai": "chai1",
    "chai1": "chai1",
    "openfold3": "openfold3",
}


def normalize_modeling_backend(modeling_backend="colabfold"):
    """
    Normalize a user-facing modeling backend name to an internal identifier.
    """
    backend_key = re.sub(r"[^a-z0-9]+", "", str(modeling_backend).lower())
    if backend_key not in _BACKEND_ALIASES:
        supported_backends = ", ".join(sorted(set(_BACKEND_ALIASES.values())))
        raise ValueError(
            "Unsupported modeling backend {!r}. Supported backends are: {}.".format(
                modeling_backend, supported_backends
            )
        )
    return _BACKEND_ALIASES[backend_key]


def backend_supports_multichain_native_inputs(modeling_backend):
    """
    Return True when the backend expects separate chains instead of legacy
    colon-joined or glycine-linked sequences.
    """
    return normalize_modeling_backend(modeling_backend) in MULTICHAIN_BACKENDS


def _ensure_output_folder(folder_path):
    folder = Path(folder_path).expanduser()
    folder.mkdir(parents=True, exist_ok=True)
    return folder


def _sanitize_label(label):
    sanitized = re.sub(r"[^A-Za-z0-9_.-]+", "-", str(label)).strip("-")
    return sanitized or "chain"


def _build_chain_records(receptor_name, receptor_seq, peptide_names, peptide_seqs):
    if len(peptide_names) != len(peptide_seqs):
        raise ValueError(
            "Peptide name and sequence lists must have the same length. Got {} names and {} sequences.".format(
                len(peptide_names), len(peptide_seqs)
            )
        )

    number_of_chains = len(peptide_names) + 1
    if number_of_chains > len(_CHAIN_IDS):
        raise ValueError(
            "APPRAISE currently supports at most {} chains per modeling job.".format(
                len(_CHAIN_IDS)
            )
        )

    chain_records = []
    for index, (chain_name, sequence) in enumerate(zip(peptide_names, peptide_seqs)):
        chain_records.append(
            {
                "chain_id": _CHAIN_IDS[index],
                "chain_name": chain_name,
                "sequence": sequence,
                "role": "peptide",
            }
        )

    chain_records.append(
        {
            "chain_id": _CHAIN_IDS[len(peptide_names)],
            "chain_name": receptor_name,
            "sequence": receptor_seq,
            "role": "receptor",
        }
    )
    return chain_records


def _build_boltz_yaml(chain_records):
    lines = ["sequences:"]
    for chain_record in chain_records:
        lines.extend(
            [
                "  - protein:",
                "      id: {}".format(chain_record["chain_id"]),
                "      sequence: {}".format(chain_record["sequence"]),
            ]
        )
    return "\n".join(lines) + "\n"


def _build_chai_fasta(chain_records):
    lines = []
    for chain_record in chain_records:
        chain_name = _sanitize_label(
            "{}_{}".format(chain_record["chain_id"], chain_record["chain_name"])
        )
        lines.append(">protein|name={}".format(chain_name))
        lines.append(chain_record["sequence"])
    return "\n".join(lines) + "\n"


def _build_openfold3_query(jobname, chain_records):
    query = {"queries": {jobname: {"chains": []}}}
    for chain_record in chain_records:
        chain_query = {
            "molecule_type": "protein",
            "chain_ids": chain_record["chain_id"],
            "sequence": chain_record["sequence"],
            "description": chain_record["chain_name"],
            "use_msas": True,
            "use_main_msas": True,
        }
        if len(chain_records) > 1:
            chain_query["use_paired_msas"] = True
        query["queries"][jobname]["chains"].append(chain_query)
    return query


def save_modeling_input(
    jobname,
    query_sequence,
    receptor_name,
    receptor_seq,
    peptide_names,
    peptide_seqs,
    folder_path="./input_fasta/",
    modeling_backend="colabfold",
):
    """
    Save an APPRAISE modeling input file for a specific backend.

    The APPRAISE job name is kept as the file stem so that downstream structure
    outputs can be matched back to the corresponding competition.
    """
    backend = normalize_modeling_backend(modeling_backend)
    output_folder = _ensure_output_folder(folder_path)

    if backend == "legacy_fasta":
        output_path = output_folder / "{}.fasta".format(jobname)
        output_path.write_text("> {}\n{}".format(jobname, query_sequence))
    else:
        chain_records = _build_chain_records(
            receptor_name, receptor_seq, peptide_names, peptide_seqs
        )
        if backend == "boltz1":
            output_path = output_folder / "{}.yaml".format(jobname)
            output_path.write_text(_build_boltz_yaml(chain_records))
        elif backend == "chai1":
            output_path = output_folder / "{}.fasta".format(jobname)
            output_path.write_text(_build_chai_fasta(chain_records))
        elif backend == "openfold3":
            output_path = output_folder / "{}.json".format(jobname)
            output_path.write_text(
                json.dumps(_build_openfold3_query(jobname, chain_records), indent=2)
            )
        else:
            raise ValueError("Unsupported backend {!r}".format(backend))

    print("$ Generated {}".format(output_path))
    return str(output_path)


def _expand_results_roots(results_path):
    normalized_path = os.path.expanduser(str(results_path))
    if any(character in normalized_path for character in "*?[]"):
        matches = sorted(glob.glob(normalized_path))
    else:
        matches = [normalized_path]
    return [Path(match) for match in matches]


def _iter_structure_files(root_path):
    if root_path.is_file():
        if root_path.suffix.lower() in STRUCTURE_FILE_EXTENSIONS:
            yield root_path
        return

    if not root_path.exists():
        return

    for extension in STRUCTURE_FILE_EXTENSIONS:
        for candidate in sorted(root_path.rglob("*{}".format(extension))):
            relative_parts = candidate.relative_to(root_path).parts[:-1]
            if any(
                part.startswith(_AUXILIARY_DIRECTORY_PREFIXES)
                or part.endswith(_AUXILIARY_DIRECTORY_SUFFIXES)
                for part in relative_parts
            ):
                continue
            yield candidate


def discover_structure_files(results_path, use_relaxed="auto"):
    """
    Discover APPRAISE-compatible structure files across flat or nested output
    directories from multiple modeling backends.

    Legacy AlphaFold/ColabFold outputs keep their relaxed/unrelaxed preference,
    while newer backends (Boltz-1, Chai-1, OpenFold3) are passed through as-is.
    """
    structure_files = []
    for root_path in _expand_results_roots(results_path):
        structure_files.extend(_iter_structure_files(root_path))

    structure_files = sorted({path.resolve() for path in structure_files})
    if not structure_files:
        return []

    relaxed_files = [path for path in structure_files if "_relaxed_" in path.name]
    unrelaxed_files = [path for path in structure_files if "_unrelaxed_" in path.name]
    backend_native_files = [
        path
        for path in structure_files
        if "_relaxed_" not in path.name and "_unrelaxed_" not in path.name
    ]

    if use_relaxed == "auto":
        prefer_relaxed = len(relaxed_files) > 0
    else:
        prefer_relaxed = bool(use_relaxed)

    selected_files = list(backend_native_files)
    if prefer_relaxed:
        selected_files.extend(relaxed_files or unrelaxed_files)
    else:
        selected_files.extend(unrelaxed_files or relaxed_files)

    return sorted({str(path) for path in selected_files})


def strip_backend_suffix(candidate_name):
    """
    Remove common modeling-tool suffixes from a file stem or folder name.
    """
    if candidate_name is None:
        return ""

    candidate = str(candidate_name)
    for extension in STRUCTURE_FILE_EXTENSIONS:
        if candidate.lower().endswith(extension):
            candidate = candidate[: -len(extension)]
            break

    if "_relaxed_" in candidate:
        return candidate.split("_relaxed_", 1)[0]
    if "_unrelaxed_" in candidate:
        return candidate.split("_unrelaxed_", 1)[0]

    previous = None
    while candidate != previous:
        previous = candidate
        candidate = re.sub(r"[_-](?:model|seed|sample|rank)[._-]?[A-Za-z0-9]+$", "", candidate)
        candidate = re.sub(r"[_-](?:scores|confidence|pae|pde|plddt)$", "", candidate)
        candidate = re.sub(r"^(?:pred|scores)[._-]model[._-]idx[._-]?[A-Za-z0-9]+$", "", candidate)
        candidate = re.sub(r"^model[._-]?[A-Za-z0-9]+$", "", candidate)

    return candidate


def is_appraise_job_name(candidate_name):
    """
    Return True when the candidate name looks like an APPRAISE job identifier.
    """
    candidate = strip_backend_suffix(candidate_name)
    return bool(
        candidate
        and (
            "_and_" in candidate
            or candidate.endswith("_receptor_model")
            or "_targeting_peptide_" in candidate
        )
    )


def extract_appraise_job_name(source_name):
    """
    Extract an APPRAISE job name from a model name, folder name, or structure path.
    """
    path = Path(str(source_name))
    candidates = [
        strip_backend_suffix(path.name),
        strip_backend_suffix(path.stem),
        strip_backend_suffix(path.parent.name),
        strip_backend_suffix(path.parent.parent.name),
    ]

    for candidate in candidates:
        if is_appraise_job_name(candidate):
            return candidate

    raise ValueError(
        "Could not infer an APPRAISE job name from {!r}. Keep the APPRAISE job name as the input file stem when running external structure-modeling tools.".format(
            str(source_name)
        )
    )


def parse_appraise_job_name(source_name):
    """
    Parse receptor and peptide names from an APPRAISE job name or structure path.
    """
    job_name = extract_appraise_job_name(source_name)

    if job_name.endswith("_receptor_model"):
        receptor_name = job_name[: -len("_receptor_model")]
        return receptor_name, []

    if "_targeting_peptide_" in job_name:
        receptor_name, peptide_name = job_name.split("_targeting_peptide_", 1)
        return receptor_name, [peptide_name]

    receptor_name, peptide_blob = job_name.split("_and_", 1)
    list_peptide_name = peptide_blob.split("_vs_")
    return receptor_name, list_peptide_name
