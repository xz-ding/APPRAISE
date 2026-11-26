# Run the β₂AR nanobody example

Choose **Runtime → Run all** in the candidate notebook. The default reads ten archived AlphaFold2-multimer structures, checks their three chain sequences, measures the two candidates in each structure and displays a rank table and heatmap. Nb6B9 ranks first and Nb80 second. The download buttons and the printed output directory provide the results. [Open the notebook](https://colab.research.google.com/github/xz-ding/APPRAISE/blob/c2c429fa2de7be3458d203a465da7e14436b5896/notebooks/Colab_APPRAISE_candidate.ipynb)

## Use the same workflow locally

Install this candidate in a Python 3.12 virtual environment. Commands write one JSON result to stdout and progress to stderr.

```sh
python -m pip install 'appraise[structures] @ git+https://github.com/xz-ding/APPRAISE.git@9c7e2e9590d90b02c46d2cf00e3be904867452ee'
python -m appraise score --example b2ar-saved --output runs/saved
python -m appraise resume --run-dir runs/saved
python -m appraise score --example b2ar-reference --output runs/reference
```

The equivalent Python interface can inspect the workload before running.

```python
from appraise.runner import RunConfig, preflight, run

config = RunConfig(example="b2ar-saved", output_dir="runs/saved")
print(preflight(config)["expected_predictions"])
result = run(config)
print(result["status"], result.get("artifacts", {}))
```

The historical reference uses 36 ordered inputs, five archived structures per input, and two observations per structure. Its Spearman correlation against the experimental order is 31/35. The two-candidate demonstration has two reciprocal inputs, five structures per input and ten independent structures per directed matchup. Gray heatmap cells mark unscheduled comparisons. Win fractions use the observed opponents. Tied final scores receive average ranks. A ranking expresses relative binding propensity for this panel. [Experiment](https://pmc.ncbi.nlm.nih.gov/articles/PMC11184338/#fig3) · [Original data](https://data.caltech.edu/records/kxjgj-tfk18)

## Select a subset from a candidate CSV

The new runner accepts the manuscript's `peptide_name` and `peptide_seq` columns for subsets of the six published nanobodies. For example, create a two-row CSV in the same format and preview both reciprocal jobs.

```python
import csv
from appraise.nanobody import panel
from appraise.runner import RunConfig, preflight

sequences = panel()["candidates"]
with open("my-nanobodies.csv", "w", newline="") as stream:
    writer = csv.writer(stream)
    writer.writerow(["peptide_name", "peptide_seq"])
    for name in ("Nb6A9", "Nb6A5"):
        writer.writerow([name, sequences[name]])

config = RunConfig(example="b2ar-panel", engine="esmfold2-fast",
                   candidates_csv="my-nanobodies.csv", output_dir="runs/my-panel")
print(preflight(config)["expected_predictions"])  # 10
```

In an environment prepared with that engine, the matching command is

```sh
python -m appraise run --example b2ar-panel --candidates-csv my-nanobodies.csv --engine esmfold2-fast --output runs/my-panel
```

The preset checks the 314-residue activated β₂AR construct and both 120-residue nanobodies. It uses the N terminus as the membrane anchor. Dmax is the maximum dimension of the receptor in ångströms. The axial ratio describes the major and minor axes of its fitted ellipsoid. The published values are 93.8 Å and 2.23, with a reported minor radius of 21 Å. The two-candidate structures have 554 residues and three protein chains. Other receptor constructs require their own geometry and justified anchor. [Construct and measurement sources](../appraise/examples/b2ar/SOURCE.md)

## Fresh predictions and interpretation

Select a free GPU and the desired validation workflow in the notebook. The two-nanobody run requests ten new structures. Six nanobodies request 180, including both copies in each self-pair. The notebook displays these counts before installing a modeling engine. Enable `save_to_drive` to keep checkpoints in `MyDrive/APPRAISE/b2ar` across Colab sessions. With that option cleared, download the results before disconnecting. Choose the same workflow and panel after reconnecting to resume a saved run. Select a different output directory for a different engine, input or scientific protocol.

| Route | Current evidence |
|---|---|
| Saved AlphaFold2 structures | The installed wheel reproduces the two-candidate ranking on CPU. |
| Archived six-nanobody measurements | The installed wheel reproduces the paper's AlphaFold2 order and correlation. |
| ColabFold2 preview ESMFold2-Fast | The v3.1.20 port with a streamed ESM-C loader ran both orders of the 554-residue competition on a free T4. All ten structures and 20 measurements passed. The five-output protocol uses seeds 0–4 with one sample per seed. [Live validation](implementation-status.md#free-colab-validation) |
| ColabFold 1.6.3 | AlphaFold2-multimer v2 with templates completed both input orders on a free T4, with ten distinct structures and 20 measurements. The five-output protocol uses five model checkpoints at one seed. Resume reused the first saved structure. [Live validation](implementation-status.md#free-colab-validation) |
| Legacy AlphaFold2 choices, ESMFold v1, Boltz-1, Chai-1 and external OpenFold3 | The separate manual notebook retains these paths. Boltz-1, Chai-1 and OpenFold3 retain their existing beta status. |

New folding with `auto` reports a blocked status while hardware qualification is pending. Saved structures and archived measurements are available on CPU. The six-candidate ESMFold2 comparison requires all 180 structures, 360 measurements, Spearman ρ of at least 29/35 and unique first and last places for Nb6B9 and Nb80. [ColabFold release](https://github.com/sokrypton/ColabFold/releases/tag/v1.6.3) · [ColabFold2 preview](https://github.com/sokrypton/ColabFold/blob/eb061aa6706959fdffe15e86d27ca08083d6ed8d/ColabFold2_preview.ipynb) · [ESMFold2-Fast model card](https://huggingface.co/biohub/ESMFold2-Fast)

## Files and failures

Each run writes `run_manifest.json`, `result.json` and `run.log`. Completed rankings include `measurements.csv`, `ranking.csv` and `heatmap.png` with SHA-256 hashes. Per-model checkpoint files retain the source identity, observed chain IDs, checksum and both candidate measurements. Resuming checks the source and measurement bytes before reusing them. Individual prediction logs and setup logs identify installation, engine or quantification failures.

The command exits with 0 for a complete run, 1 for failure, 3 for an external or qualification block and 130 for an interruption. An incomplete result identifies its stage, active input, counts, error and next action. Invoke `resume` after restoring the dependency or destination. A changed configuration gets a separate run directory. Historical measurements include only the first 70 residues in their sequence field. Newly quantified structures record all 120.
