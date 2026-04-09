import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from appraise.input_fasta_prep import get_complex_fastas
from appraise.modeling_backends import (
    discover_structure_files,
    extract_appraise_job_name,
    normalize_modeling_backend,
    parse_appraise_job_name,
)
from appraise.pymol_quantify_peptide_binding import (
    generate_pdb_path_list,
    parse_pdb_file_name,
)
from appraise.utilities import get_peptide_list_from_model_names


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

    def test_job_name_parsing_supports_legacy_and_new_backends(self):
        legacy_name = "LY6A_and_pepA_vs_pepB_unrelaxed_rank_001_model_1_seed_000.pdb"
        boltz_name = "/tmp/predictions/LY6A_and_pepA_vs_pepB/LY6A_and_pepA_vs_pepB_model_0.cif"
        chai_name = "/tmp/chai/LY6A_and_pepA_vs_pepB/pred.model_idx_0.cif"

        for source_name in [legacy_name, boltz_name, chai_name]:
            self.assertEqual(extract_appraise_job_name(source_name), "LY6A_and_pepA_vs_pepB")
            self.assertEqual(parse_appraise_job_name(source_name), ("LY6A", ["pepA", "pepB"]))
            self.assertEqual(parse_pdb_file_name(source_name), ("LY6A", ["pepA", "pepB"]))


class ExistingPipelineCompatibilityTests(unittest.TestCase):
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
