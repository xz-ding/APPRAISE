import json
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK_PATH = REPO_ROOT / "Colab_APPRAISE.ipynb"


class ColabNotebookRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.notebook = json.loads(NOTEBOOK_PATH.read_text())
        cls.cell_sources = ["".join(cell.get("source", [])) for cell in cls.notebook["cells"]]

    def _find_cell(self, expected_text):
        for source in self.cell_sources:
            if expected_text in source:
                return source
        self.fail(f"Could not find notebook cell containing: {expected_text!r}")

    def test_step0_installs_appraise_from_pypi(self):
        self._find_cell("### **0.2 Install APPRAISE package**")
        install_cell = self._find_cell("!pip install -q --upgrade appraise")

        self.assertIn("Install the latest released APPRAISE package from PyPI", install_cell)
        self.assertIn("recommended path for the public notebook", install_cell)
        self.assertIn("replace this command with a GitHub install manually", install_cell)
        self.assertNotIn("git+https://github.com/xz-ding/APPRAISE.git", install_cell)

    def test_step1_backend_selector_lists_new_backends(self):
        source = self._find_cell("modeling_backend = 'colabfold'")
        self.assertIn("['colabfold', 'boltz1', 'chai1', 'openfold3']", source)
        self.assertIn("Boltz-1 uses YAML inputs", source)
        self.assertIn("Chai-1 uses multichain FASTA inputs", source)
        self.assertIn("OpenFold3 uses JSON query files", source)

    def test_step1_passes_backend_into_input_preparation(self):
        source = self._find_cell("get_complex_fastas(")
        self.assertIn("modeling_backend=modeling_backend", source)

    def test_step2_includes_boltz_and_chai_blocks(self):
        boltz_step = self._find_cell("### **Step 2C - Predict structures with Boltz-1**")
        chai_step = self._find_cell("### **Step 2D - Predict structures with Chai-1**")
        boltz_run = self._find_cell("boltz_executable = shutil.which('boltz')")
        chai_run = self._find_cell("chai_executable = shutil.which('chai-lab')")

        self.assertIn("official `boltz predict` CLI", boltz_step)
        self.assertIn("**BETA support:**", boltz_step)
        self.assertIn("official `chai-lab fold` CLI", chai_step)
        self.assertIn("**BETA support:**", chai_step)
        self.assertIn("glob.glob(os.path.join(input_dir, '*.yaml'))", boltz_run)
        self.assertIn("command = [", boltz_run)
        self.assertIn("'predict'", boltz_run)
        self.assertIn("glob.glob(os.path.join(input_dir, '*.fasta'))", chai_run)
        self.assertIn("command = [chai_executable, 'fold']", chai_run)

    def test_step2a_preserves_legacy_alphafold_multimer_flow(self):
        settings = self._find_cell('model_type = "alphafold2_multimer_v2"')
        run_cell = self._find_cell("download_alphafold_params(selected_model_type, Path(default_data_dir))")

        self.assertIn(
            '["alphafold2_multimer_v1", "alphafold2_multimer_v2", "alphafold2_multimer_v3"]',
            settings,
        )
        self.assertIn("keep `modeling_backend='colabfold'` in Step 1.4", self._find_cell("### **Step 2A - Predict structures with AlphaFold-multimer**"))
        self.assertIn("from colabfold.batch import get_queries, run", run_cell)
        self.assertIn("num_relax = 1 if use_amber else 0", run_cell)
        self.assertIn('user_agent="colabfold/appraise-colab"', run_cell)

    def test_step2b_preserves_legacy_esmfold_flow(self):
        step = self._find_cell("### **Step 2B - Predict structures with ESMFold**")
        install = self._find_cell("Transformers-based ESMFold workflow")
        run_cell = self._find_cell("pdb_filename = result_dir.joinpath(f\"{jobname}_unrelaxed_ptm")

        self.assertIn("keep `modeling_backend='colabfold'` in Step 1.4", step)
        self.assertIn("def safe_filename", install)
        self.assertIn("def setup_logging", install)
        self.assertIn("default_data_dir = Path.home() / '.cache' / 'huggingface'", install)
        self.assertIn('tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")', self._find_cell('tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")'))
        self.assertIn("if 'logging_setup' not in globals():", run_cell)
        self.assertIn('with open(pdb_filename,"w") as out:', run_cell)

    def test_step2_guide_mentions_external_openfold3(self):
        source = self._find_cell("**Guide to step 2**")
        self.assertIn("`openfold3` for external OpenFold3 runs followed by Step 3", source)
        self.assertIn("do **not** embed OpenFold3 directly", source)
        self.assertIn("currently **BETA** only", source)

    def test_step3_documents_modern_structure_outputs(self):
        source = self._find_cell("**Guide to Step 3**")
        self.assertIn("nested `.cif` and `.mmcif` outputs", source)
        self.assertIn("currently **BETA**", source)
        notes = self._find_cell("<font color=\"grey\">**Notes:**")
        self.assertIn("Support beyond AlphaFold-multimer and ESMFold is currently **BETA**", notes)
        self.assertIn("usually do **not** need to rename structure outputs", notes)
        self.assertNotIn("14 characters", notes)

    def test_step32_surfaces_structure_discovery_and_script_output(self):
        step32 = self._find_cell("APPRAISE> Found {len(structure_files)} structure files")
        self.assertIn("from appraise.modeling_backends import discover_structure_files", step32)
        self.assertIn("No APPRAISE-compatible structure files were found", step32)
        self.assertIn("except subprocess.CalledProcessError as exc", step32)
        self.assertIn("APPRAISE Step 3.2 failed", step32)
        self.assertIn("capture_output=True", step32)
        self.assertIn("APPRAISE> Step 3.2 completed successfully.", step32)

    def test_step4_reports_empty_database_more_helpfully(self):
        pairwise_cell = self._find_cell("No rows in the selected database match receptor_of_interest")
        self.assertIn("read_appraise_database(database_path)", pairwise_cell)
        self.assertIn("0 peptide variants were detected after parsing model names", pairwise_cell)
        self.assertIn("contains only receptor-only models", pairwise_cell)


if __name__ == "__main__":
    unittest.main()
