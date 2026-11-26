import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd
import pytest

from appraise.input_fasta_prep import get_complex_fastas
from appraise.legacy_esmfold import restore_esmfold_v1_chains
from appraise.modeling_backends import (
    discover_structure_files,
    extract_appraise_job_name,
    normalize_modeling_backend,
    parse_appraise_job_name,
)
from appraise.pymol_quantify_peptide_binding import (
    count_database_measurements,
    generate_pdb_path_list,
    parse_pdb_file_name,
)
from appraise.utilities import get_peptide_list_from_model_names


def test_esmfold_v1_recovery_preserves_natural_glycines_and_coordinates():
    first = "AG" + "G" * 24 + "CG"
    second = "GCGA"
    joined = first + "G" * 30 + second
    amino_acid = {"A": "ALA", "C": "CYS", "G": "GLY"}
    original = [f"ATOM  {i:5d}  CA  {amino_acid[aa]} A{i:4d}    "
                f"{i:8.3f}{0:8.3f}{0:8.3f}{1:6.2f}{89:6.2f}           C"
                for i, aa in enumerate(joined, 1)]
    pdb = "\n".join(original + ["TER", "END", ""])
    repaired = [line for line in restore_esmfold_v1_chains(pdb, [first, second]).splitlines()
                if line.startswith("ATOM")]
    retained = original[:len(first)] + original[-len(second):]
    assert [line[30:] for line in repaired] == [line[30:] for line in retained]
    assert [line[17:20] for line in repaired] == [line[17:20] for line in retained]
    assert [line[21] for line in repaired] == ["A"] * len(first) + ["B"] * len(second)
    assert [int(line[22:26]) for line in repaired if line[21] == "B"] == [1, 2, 3, 4]
    with pytest.raises(ValueError, match="original separated"):
        restore_esmfold_v1_chains(pdb, None)
    with pytest.raises(ValueError, match="does not match"):
        restore_esmfold_v1_chains(pdb, [second, first])


class ModelingInputCompatibilityTests(unittest.TestCase):
    def test_legacy_backend_still_writes_colon_joined_fasta(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            query_sequences, jobnames = get_complex_fastas(
                receptor_name="LY6A",
                receptor_seq="RRR",
                list_peptide1_names=["pepA"],
                list_peptide1_seqs=["AAA"],
                mode="single",
                folder_path=tmpdir,
                modeling_backend="colabfold",
            )

            self.assertEqual(query_sequences, ["AAA:RRR", "RRR"])
            self.assertEqual(jobnames, ["LY6A_and_pepA", "LY6A_receptor_model"])

            fasta_path = Path(tmpdir) / "LY6A_and_pepA.fasta"
            self.assertTrue(fasta_path.exists())
            self.assertEqual(fasta_path.read_text(), "> LY6A_and_pepA\nAAA:RRR")

    def test_versioned_legacy_backend_aliases_still_write_legacy_fastas(self):
        legacy_aliases = [
            "alphafold2_multimer",
            "alphafold2_multimer_v2",
            "alphafold2_multimer_v3",
            "esmfold",
            "esmfold2",
        ]

        for alias in legacy_aliases:
            with self.subTest(alias=alias):
                with tempfile.TemporaryDirectory() as tmpdir:
                    query_sequences, _ = get_complex_fastas(
                        receptor_name="LY6A",
                        receptor_seq="RRR",
                        list_peptide1_names=["pepA"],
                        list_peptide1_seqs=["AAA"],
                        mode="single",
                        folder_path=tmpdir,
                        modeling_backend=alias,
                        use_glycine_linker=True,
                        glycine_linker_length=5,
                    )

                    self.assertEqual(normalize_modeling_backend(alias), "legacy_fasta")
                    self.assertEqual(query_sequences, ["AAAGGGGGRRR", "RRR"])

                    fasta_path = Path(tmpdir) / "LY6A_and_pepA.fasta"
                    self.assertEqual(
                        fasta_path.read_text(),
                        "> LY6A_and_pepA\nAAAGGGGGRRR",
                    )
                    records = json.loads(fasta_path.with_suffix(".chains.json").read_text())
                    self.assertEqual(
                        [(r["chain_id"], r["role"], r["sequence"]) for r in records],
                        [("A", "peptide", "AAA"), ("B", "receptor", "RRR")],
                    )

    def test_boltz_input_generation_uses_multichain_yaml(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            query_sequences, _ = get_complex_fastas(
                receptor_name="LY6A",
                receptor_seq="RRR",
                list_peptide1_names=["pepA"],
                list_peptide1_seqs=["AAA"],
                mode="single",
                folder_path=tmpdir,
                modeling_backend="boltz1",
                use_glycine_linker=True,
            )

            self.assertEqual(query_sequences, ["AAA:RRR", "RRR"])

            fasta_path = Path(tmpdir) / "LY6A_and_pepA.yaml"
            self.assertEqual(
                fasta_path.read_text(),
                "sequences:\n  - protein:\n      id: A\n      sequence: AAA\n  - protein:\n      id: B\n      sequence: RRR\n",
            )

    def test_chai_input_generation_uses_multichain_fasta(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            get_complex_fastas(
                receptor_name="LY6A",
                receptor_seq="RRR",
                list_peptide1_names=["pepA"],
                list_peptide1_seqs=["AAA"],
                mode="single",
                folder_path=tmpdir,
                modeling_backend="chai1",
            )

            fasta_path = Path(tmpdir) / "LY6A_and_pepA.fasta"
            self.assertEqual(
                fasta_path.read_text(),
                ">protein|name=A_pepA\nAAA\n>protein|name=B_LY6A\nRRR\n",
            )

    def test_openfold3_input_generation_writes_query_json(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            get_complex_fastas(
                receptor_name="LY6A",
                receptor_seq="RRR",
                list_peptide1_names=["pepA"],
                list_peptide1_seqs=["AAA"],
                mode="single",
                folder_path=tmpdir,
                modeling_backend="openfold3",
            )

            json_path = Path(tmpdir) / "LY6A_and_pepA.json"
            payload = json.loads(json_path.read_text())
            self.assertIn("queries", payload)
            query = payload["queries"]["LY6A_and_pepA"]
            self.assertEqual(len(query["chains"]), 2)
            self.assertEqual(query["chains"][0]["chain_ids"], "A")
            self.assertEqual(query["chains"][0]["sequence"], "AAA")
            self.assertEqual(query["chains"][1]["chain_ids"], "B")
            self.assertEqual(query["chains"][1]["sequence"], "RRR")


class StructureDiscoveryCompatibilityTests(unittest.TestCase):
    def test_legacy_relaxed_preference_is_preserved(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            results_dir = Path(tmpdir)
            relaxed = results_dir / "LY6A_and_pepA_relaxed_rank_001_model_1_seed_000.pdb"
            unrelaxed = results_dir / "LY6A_and_pepA_unrelaxed_rank_001_model_1_seed_000.pdb"
            relaxed.write_text("RELAXED")
            unrelaxed.write_text("UNRELAXED")

            self.assertEqual(discover_structure_files(tmpdir, use_relaxed="auto"), [str(relaxed.resolve())])
            self.assertEqual(discover_structure_files(tmpdir, use_relaxed=False), [str(unrelaxed.resolve())])
            self.assertEqual(generate_pdb_path_list(tmpdir, use_relaxed="auto"), [str(relaxed.resolve())])

    def test_nested_backend_outputs_are_discovered(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            boltz_path = root / "boltz_results" / "predictions" / "LY6A_and_pepA_vs_pepB" / "LY6A_and_pepA_vs_pepB_model_0.cif"
            chai_path = root / "chai_results" / "LY6A_and_pepA_vs_pepB" / "pred.model_idx_0.cif"
            openfold_path = root / "openfold3_results" / "LY6A_and_pepA_vs_pepB_model_0.cif"
            for structure_path in [boltz_path, chai_path, openfold_path]:
                structure_path.parent.mkdir(parents=True, exist_ok=True)
                structure_path.write_text("data_mock\n#\n")

            discovered = discover_structure_files(tmpdir, use_relaxed="auto")
            expected = sorted(
                [
                    str(boltz_path.resolve()),
                    str(chai_path.resolve()),
                    str(openfold_path.resolve()),
                ]
            )
            self.assertEqual(discovered, expected)
            self.assertEqual(generate_pdb_path_list(tmpdir, use_relaxed="auto"), expected)

    def test_colabfold_auxiliary_template_directories_are_ignored(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            predicted_pdb = root / "LY6A_and_pepA_vs_pepB_unrelaxed_rank_001_model_1_seed_000.pdb"
            env_template = root / "LY6A_and_pepA_vs_pepB_env" / "templates_101" / "1abc.cif"
            pairgreedy_template = root / "LY6A_and_pepA_vs_pepB_pairgreedy" / "candidate.cif"

            predicted_pdb.write_text("MODEL     1\n")
            env_template.parent.mkdir(parents=True, exist_ok=True)
            env_template.write_text("data_mock\n#\n")
            pairgreedy_template.parent.mkdir(parents=True, exist_ok=True)
            pairgreedy_template.write_text("data_mock\n#\n")

            discovered = discover_structure_files(tmpdir, use_relaxed="auto")

            self.assertEqual(discovered, [str(predicted_pdb.resolve())])

    def test_job_name_parsing_supports_legacy_and_new_backends(self):
        legacy_name = "LY6A_and_pepA_vs_pepB_unrelaxed_rank_001_model_1_seed_000.pdb"
        boltz_name = "/tmp/predictions/LY6A_and_pepA_vs_pepB/LY6A_and_pepA_vs_pepB_model_0.cif"
        chai_name = "/tmp/chai/LY6A_and_pepA_vs_pepB/pred.model_idx_0.cif"

        for source_name in [legacy_name, boltz_name, chai_name]:
            self.assertEqual(extract_appraise_job_name(source_name), "LY6A_and_pepA_vs_pepB")
            self.assertEqual(parse_appraise_job_name(source_name), ("LY6A", ["pepA", "pepB"]))
            self.assertEqual(parse_pdb_file_name(source_name), ("LY6A", ["pepA", "pepB"]))


class ExistingPipelineCompatibilityTests(unittest.TestCase):
    def test_count_database_measurements_ignores_header_only_csv(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            database_path = Path(tmpdir) / "database_APPRAISE_measurements.csv"
            database_path.write_text("model_name,receptor_name,peptide_name\n")
            self.assertEqual(count_database_measurements(database_path), 0)

            database_path.write_text(
                "model_name,receptor_name,peptide_name\n"
                "LY6A_and_pepA_unrelaxed_rank_001,LY6A,pepA\n"
                "LY6A_and_pepB_unrelaxed_rank_001,LY6A,pepB\n"
            )
            self.assertEqual(count_database_measurements(database_path), 2)

    def test_get_peptide_list_from_model_names_handles_modern_model_names(self):
        df = pd.DataFrame(
            {
                "model_name": [
                    "LY6A_and_pepA_unrelaxed_rank_001_model_1_seed_000",
                    "LY6A_and_pepA_unrelaxed_rank_002_model_1_seed_001",
                    "pred.model_idx_0",
                    "LY6A_and_pepB_model_0",
                ],
                "peptide_name": ["pepA", "pepA", "pepB", "pepB"],
                "peptide_seq": ["AAA", "AAA", "BBB", "BBB"],
            }
        )
        df.loc[2, "model_name"] = "LY6A_and_pepB_model_0"

        peptide_names, peptide_seqs = get_peptide_list_from_model_names(df)

        self.assertEqual(set(peptide_names), {"pepA", "pepB"})
        self.assertEqual(set(peptide_seqs), {"AAA", "BBB"})


if __name__ == "__main__":
    unittest.main()
