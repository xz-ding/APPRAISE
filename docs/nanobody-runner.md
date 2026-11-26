# Run the β₂AR nanobody example

Choose **Runtime → Run all** in the candidate notebook. The default reads ten archived AlphaFold2-multimer structures, checks their three chain sequences, measures the two candidates in each structure and displays a rank table and heatmap. Nb6B9 ranks first and Nb80 second. The download buttons and the printed output directory provide the results. [Open the notebook](https://colab.research.google.com/github/xz-ding/APPRAISE/blob/31934881b32c2ba9c101084bcffad972c0e0f3e6/notebooks/Colab_APPRAISE_candidate.ipynb)

## Use the same workflow locally

Install this candidate in a Python 3.12 virtual environment. Commands write one JSON result to stdout and progress to stderr.

```sh
python -m pip install 'appraise[structures] @ git+https://github.com/xz-ding/APPRAISE.git@1f6de6b18250ff387bf261a72eec65121a551027'
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
| JapanFold API folding | Planned for the Colab notebook and Python interface. [Integration requirements](#planned-japanfold-integration) |

### Recover joined ESMFold v1 structures

For an older single-chain ESMFold v1 PDB, use the original separated sequences in the same order as the model inputs. Newly generated glycine-joined FASTA inputs carry that information in a `.chains.json` file. Recover the chain boundaries before running quantification, keeping the original PDB and writing the recovered file to a separate results folder.

```python
import json
from pathlib import Path
from appraise.legacy_esmfold import restore_esmfold_v1_chains

source = Path("old-results/example_and_candidate.pdb")
original = json.loads(Path("inputs/example_and_candidate.chains.json").read_text())
restored = restore_esmfold_v1_chains(source.read_text(), [chain["sequence"] for chain in original])
destination = Path("recovered-results") / source.name
destination.parent.mkdir(exist_ok=True)
destination.write_text(restored)
```

When no sidecar exists, supply the original separated sequence strings directly. The recovery step checks the exact joined sequence and removes the intervening artificial glycines. If the original sequences are unavailable, the chain assignment is ambiguous.

### Estimated ESMFold2 usage costs

| Where prediction runs | Estimated charge and who bills it |
|---|---|
| ColabFold2 ESMFold2-Fast on free Colab | **$0 estimated GPU charge** while Google assigns a free runtime. Google sets Colab's GPU usage limits and charges for optional paid Colab compute. [Google Colab FAQ](https://research.google.com/colaboratory/faq.html#usage-limits) |
| JapanFold ESMFold2-Fast | **Estimated $0.0011 per prediction or $0.20 for 180 similar predictions.** The estimate uses the $0.194861 JapanFold charged for 180 predictions of 554-residue complexes. JapanFold, operated by ai&, bills actual occupied processor time to the JapanFold account. Its rate was $0.26 per processor-hour on November 26, 2026. [JapanFold pricing and limits](https://japanfold.aiand.com/docs/models-and-limits/) |

The first hosted prediction cost $0.002617. Multiply estimated processor-hours by the provider's live rate when planning a new run. Actual JapanFold charges vary by input, batching, model settings, and the processor time shown on each job's usage receipt. These estimates apply to the Fast checkpoint and the benchmark settings described in the validation record. [APPRAISE validation record](implementation-status.md)

For fresh β₂AR predictions, `auto` selects AlphaFold2-multimer v2 with templates on a Colab Tesla T4 with at least 4 GiB of free GPU memory. The threshold leaves headroom above the 3,083 MiB observed in the guided run. A resumed run keeps the selected engine. Saved structures and archival measurements are available on CPU. ESMFold2-Fast remains experimental. The all-JapanFold comparison reached ρ = 33/35 against the experimental order and placed Nb80 fifth. The six-candidate qualification gate requires 180 structures, 360 measurements, ρ of at least 29/35, and unique first and last places for Nb6B9 and Nb80. [Benchmark results](implementation-status.md#japanfold-benchmark-and-ranking)

## Planned JapanFold integration

Use the same APPRAISE workflow from Colab and Python. The user supplies a JapanFold API key through Colab Secrets or the `JAPANFOLD_API_KEY` environment variable, a maximum number of parallel jobs, and a spending cap in USD. Before submission, show the folding method, prediction count, and labeled estimate of JapanFold's charge alongside the spending cap. Each target counts toward the parallel limit, including targets grouped into one provider job. Reserve the possible charge for each outstanding prediction and reconcile it against the provider receipt before starting more work. [JapanFold prediction API](https://japanfold.aiand.com/docs/predictions/) · [JapanFold billing](https://japanfold.aiand.com/docs/billing/)

Save submission identities and JapanFold job IDs so a resumed run monitors existing jobs without duplicate charges. Show progress, partial failures, and recorded spend in APPRAISE. Download results automatically and verify each structure's chains, sequence, and sample identity before passing it to quantification and scoring. Start the final ranking when the selected comparison has the required coverage. Verify each supported folding method's actual settings and output format before listing it as available. [JapanFold jobs API](https://japanfold.aiand.com/docs/jobs/)

## Files and failures

Each run writes `run_manifest.json`, `result.json` and `run.log`. Completed rankings include `measurements.csv`, `ranking.csv` and `heatmap.png` with SHA-256 hashes. Per-model checkpoint files retain the source identity, observed chain IDs, checksum and both candidate measurements. Resuming checks the source and measurement bytes before reusing them. Individual prediction logs and setup logs identify installation, engine or quantification failures.

The default score matrix uses blue for positive scores, red for negative scores, and a neutral midpoint at zero.

The command exits with 0 for a complete run, 1 for failure, 3 for an external or qualification block and 130 for an interruption. An incomplete result identifies its stage, active input, counts, error and next action. Invoke `resume` after restoring the dependency or destination. A changed configuration gets a separate run directory. Historical measurements include only the first 70 residues in their sequence field. Newly quantified structures record all 120.
