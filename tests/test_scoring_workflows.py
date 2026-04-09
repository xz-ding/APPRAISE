import tempfile
import unittest
from pathlib import Path

import matplotlib
import pandas as pd

from appraise.input_fasta_prep import get_complex_fastas
from appraise.score_calculation import calculate_scores
from appraise.utilities import plot_heatmap, sort_df_by_peptides_and_cleanup


matplotlib.use("Agg")


REPO_ROOT = Path(__file__).resolve().parents[1]


class ScoringWorkflowTests(unittest.TestCase):
    def test_fasta_generation_count_matches_demo_expectation(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            peptide_list_path = (
                REPO_ROOT
                / "demo"
                / "example_input_sequences_from_manuscript"
                / "PD-L1_5x5"
                / "PD-L1_5x5_peptide_list_corrected.csv"
            )
            df_peptides = pd.read_csv(peptide_list_path)
            query_sequences, jobnames = get_complex_fastas(
                receptor_name="PDL1",
                receptor_seq="ACDEFGHIKLMNPQRSTVWY",
                list_peptide1_names=df_peptides["peptide_name"].tolist(),
                list_peptide1_seqs=df_peptides["peptide_seq"].tolist(),
                mode="pairwise",
                folder_path=tmpdir,
            )

            self.assertEqual(len(query_sequences), 26)
            self.assertEqual(len(jobnames), 26)
            self.assertTrue((Path(tmpdir) / "PDL1_receptor_model.fasta").exists())

    def test_pairwise_scoring_and_heatmap_still_work(self):
        database_path = (
            REPO_ROOT
            / "demo"
            / "demo_100AAV_screening"
            / "stage_2"
            / "database_APPRAISE_measurements_stage_2_example.csv"
        )
        df = pd.read_csv(database_path)
        df["receptor_Rminor"] = 46.68 / 1.74 / 2

        scored = calculate_scores(df.copy(), version=1.2)
        peptide_names = scored["peptide_name"].dropna().unique().tolist()
        sorted_df = sort_df_by_peptides_and_cleanup(scored.copy(), peptide_names)
        df_average = (
            sorted_df.groupby(["peptide_name", "competitor"]).mean(numeric_only=True).reset_index()
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            current_dir = Path.cwd()
            try:
                Path(tmpdir).mkdir(exist_ok=True)
                import os

                os.chdir(tmpdir)
                ranking, _, _ = plot_heatmap(
                    df_average.copy(),
                    feature_of_interest="Delta_B",
                    receptor_of_interest="LY6A",
                    save_figure=True,
                )
            finally:
                os.chdir(current_dir)

            self.assertEqual(ranking[:5], ["Dis90", "SRK-39", "SRK-74", "PHP.B", "SRK-50"])
            self.assertTrue((Path(tmpdir) / "LY6A_ranked_by_Delta_B.png").exists())

    def test_pooled_scoring_still_works(self):
        database_path = (
            REPO_ROOT
            / "demo"
            / "demo_100AAV_screening"
            / "stage_1"
            / "database_APPRAISE_measurements_stage_1_example.csv"
        )
        df = pd.read_csv(database_path)
        df["receptor_Rminor"] = 46.68 / 1.74 / 2

        scored = calculate_scores(df.copy(), version=1.2, pairwise=False)

        self.assertEqual(scored["peptide_name"].nunique(), 100)
        top_hits = (
            scored.groupby("peptide_name")
            .mean(numeric_only=True)
            .sort_values("B_POI", ascending=False)
            .head(3)
            .index.tolist()
        )
        self.assertEqual(top_hits, ["SRK-53", "SRK-50", "Dis90"])


if __name__ == "__main__":
    unittest.main()
