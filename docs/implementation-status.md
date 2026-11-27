# APPRAISE 1.3.2 release record

The November 25 owner-supplied plan is recorded byte for byte in `AGENT_EXECUTION_PLAN.md`. Its SHA-256 is `3337349356d7503966539223616ed40017c56cf653a7e9ae355038ca650a66e2`. Work started from `8898c203d35b51d64b4766342ed9a9a32630cc13` on canonical `main`. The package and notebook merged in separate pull requests. [Package PR](https://github.com/xz-ding/APPRAISE/pull/1) · [Notebook PR](https://github.com/xz-ding/APPRAISE/pull/2) · [Machine-readable checkpoint](implementation-state.json)

## What runs today

The notebook, noninteractive command and Python interface share the published β₂AR workflow. Each sample records the producer, input role and copy, observed chain IDs, structure and measurement hashes. Recovery tests cover interruptions, damaged files and changed chain metadata. The preview adapter requests one ESMFold2-Fast diffusion sample per seed and the fallback requests one AlphaFold2 checkpoint per prediction. The guided test case contains three protein chains and 554 residues. [Public notebook](https://colab.research.google.com/github/xz-ding/APPRAISE/blob/main/notebooks/Colab_APPRAISE.ipynb) · [Usage](nanobody-runner.md)

| Evidence | Result |
|---|---|
| Published six-nanobody measurements | All 36 ordered inputs, 180 structures and 360 measurements represented. The APPRAISE 1.2 order is Nb6B9, Nb6A5, Nb6A9, Nb6A10, Nb6B11, Nb80. Spearman ρ against the experimental order is 31/35. |
| CPU structure-to-ranking run | Headless PyMOL measured the full chains from both directions in ten archived structures. All 20 perspectives were retained. Nb6B9 ranked first and Nb80 second. |
| Package | The wheel and source archive are built from a clean source export. The release check compares their package contents with source, installs the wheel outside the checkout, runs both examples and verifies resume. The hosted job retains the exact artifacts and their hashes. |
| Live free Colab | The built-in browser connected to a free Tesla T4 with 15 GiB VRAM and 12.67 GiB system RAM. Temporary checkpoint write and read-back passed. |
| AlphaFold2-multimer v2 | ColabFold 1.6.3 used templates and completed both input orders. It produced ten distinct 554-residue structures, 20 measurements, the ranked table and heatmap. Nb6B9 placed first. |
| ESMFold2-Fast | The v3.1.20 port with the streamed ESM-C loader completed both input orders. It produced ten distinct 554-residue structures, 20 measurements, the ranked table and heatmap. Nb6B9 placed first. |
| Six-nanobody comparison | JapanFold ESMFold2-Fast completed 36 ordered inputs, 180 distinct structures and 360 measurements. Spearman ρ was 33/35 against the experimental ranking. Nb80 finished fifth. The qualification gate requires sixth place. [Ranking](results/esmfold2-japanfold-b2ar-six-ranking.csv) |

The latest candidate passed 44 scoped tests and five subtests locally. The earlier candidate passed 40 scoped tests and five subtests on GitHub’s Linux runner. They cover the two notebooks, backend I/O, the nanobody workflow, and archived measurement parsing. They also test ties, missing coverage, ensemble identities, checksum recovery, and signed colors in the guided and manual score matrices. [Hosted CPU checks](https://github.com/xz-ding/APPRAISE/actions/runs/45877547525)

## Free-Colab validation

The original port exhausted host memory while loading ESM-C. Streaming the same weight records resolved that stage. The original model graph expected 447 features from a 451-feature checkpoint. The pinned v3.1.20 port contains the upstream shape fix. The first new ESMFold2-Fast model and the first AlphaFold2-multimer model both passed written-chain and native-confidence checks. [Upstream checkpoint fix](https://github.com/sokrypton/alphafold3/commit/d87fa1eacaf275847f69e7f3ea854c362985de18)

The ESMFold2-Fast run took 58 minutes and 53 seconds after setup, and sampled peak GPU memory was 11,333 MiB. All ten structures passed chain and hash verification. The pinned environment was newly installed in the session and reused the model weights and ESM-C disk cache built earlier in that runtime. The table, measurements, heatmap and manifest hashes match the live notebook. [Ranking CSV](results/colab-t4-esmfold2-fast-b2ar-two-ranking.csv) · [Runtime and file hashes](implementation-state.json)

The AlphaFold2-multimer v2 test resumed from the first saved structure using ColabFold alignments and templates prepared for the exact two inputs. The resumed segment took 36 minutes and 18 seconds, excluding this preparation. Sampled peak GPU memory was 3,083 MiB. After Colab reconnected to the same runtime, the first structure's checksum still matched and all ten structures and final artifacts passed verification. [Ranking CSV](results/colab-t4-alphafold2-multimer-b2ar-two-ranking.csv) · [Runtime and file hashes](implementation-state.json)

The T4 handled the measured 554-residue three-chain workload. The guided ESMFold2 timing implied about 17 T4 GPU hours for the 180-prediction panel. Colab's free usage limit interrupted the run after 18 structures and 36 measurements.

The six-nanobody run used `/content/appraise-local/b2ar/b2ar-six-esmfold2-fast` on the GPU runtime. The automatic approval review declined the Drive filesystem mount because it requested account-wide file access. Six checksum-verified checkpoint archives were saved through the Drive web editor under `My Drive/APPRAISE/b2ar`. The fresh CPU runtime verified all 45 restored structure, confidence and measurement files across the 15 saved samples. Those samples form the durable T4 checkpoint.

On November 26, Colab again declined a free T4 allocation because of the account's GPU usage limit. The pinned APPRAISE `resume` command on a fresh Colab CPU runtime reused three ordered inputs, 15 structures and 30 measurements under the original run ID. All 45 file hashes and the 15 completed sample records still matched. The command exited at prediction 16 with `gpu_unavailable`.

For a T4-only replay, restore the six archives from `001-001` through `011-015` in ascending order into the same run directory. Each Markdown file starts with a marker and JSON metadata, followed by a base64-encoded gzip archive and an end marker. Verify each archive against the SHA-256 in the execution record before unpacking it. The last archive supplies the run manifest for all 15 saved samples. Install APPRAISE at `9c7e2e9590d90b02c46d2cf00e3be904867452ee` with the tested ESMFold2-Fast environment and run `python -m appraise resume --run-dir /content/appraise-local/b2ar/b2ar-six-esmfold2-fast`. The saved run ID is `1771fadf-bdf3-422a-bec7-c7c12b8980bd` and the next prediction is 16. [Archive checksums](implementation-state.json)

The guided runner covers the exact published nanobodies, their subsets and their verified receptor construct. The small-protein checks below establish runtime and format evidence for two additional inputs. Broader custom receptors remain outside the guided runner's measured scope. Synthetic tests verify legacy ESMFold v1 export and recovery from exact source sequences. The manual notebook retains the existing advanced engines. The SKEMPI 2.0 follow-up and its within-target curation rules are in section 7 of the execution plan.

## Peptide and miniprotein runtime checks

On November 27, a free Colab Tesla T4 ran the released APPRAISE 1.3.2 adapter with the pinned ColabFold2 ESMFold2-Fast port 3.1.20. The T4 reported 15,360 MiB total VRAM and 14,913 MiB free at setup. System RAM was 12.67 GiB. The runtime used Python 3.13.16, NVIDIA driver 580.82.07, and driver-reported CUDA 13.0. Setup took 63.15 seconds, including 11.69 seconds for package installation and 42.66 seconds for the pinned weight and script cache. [Runtime and dependency details](implementation-state.json)

The test harness passed two deposited folding controls and the published 314-residue β₂AR construct to the released prediction and quantification functions. Chains A and B held the controls and chain C held the receptor. Each ordered input requested one diffusion sample at seed 0, ten recycles and fifteen sampling steps, with query-only input, int8 weights and XLA attention. [Chignolin](https://www.rcsb.org/structure/1UAO) · [Trp-cage TC5b](https://www.rcsb.org/structure/1L2Y) · [Villin HP35 N68H](https://www.rcsb.org/structure/1YRF) · [Crambin](https://www.rcsb.org/structure/1CRN) · [β₂AR sequence source](../appraise/examples/b2ar/SOURCE.md)

| Case | Control lengths | Total residues | Prediction command | Logged inference | Peak VRAM at 1 s | APPRAISE rows |
|---|---:|---:|---:|---:|---:|---:|
| Peptide | 10 + 20 | 344 | 542.25 s | 84.29 s | 11,331 MiB | 2 |
| Miniprotein | 35 + 46 | 395 | 160.12 s | 126.57 s | 11,331 MiB | 2 |

The prediction command includes model startup and feature preparation. Feature preparation, including ESM-C, took 442.50 seconds for the first case and 22.08 seconds for the later case. Both completed with XLA attention after the upstream optional Tokamax autotuner reported an MLIR API error.

Released APPRAISE reproduced four measurement rows and verified the exact chain sequences, finite coordinates and native pLDDT. Both archives contained complete PAE and contact probability arrays in chain order. Their SHA-256 checksums matched the live Colab exports. The measured points for the Fast port on T4 are 344 and 395 residues here and 554 residues in the earlier nanobody example. These are runtime and output-format measurements. The recommendation continues to follow the six-nanobody ranking criterion. This follow-up cost $0, bringing tracked plan spending to $0.194861 of the $20 cap. [Inputs, artifact hashes and timings](implementation-state.json)

## JapanFold benchmark and ranking

The first comparison combined 15 ColabFold2 preview structures with 165 JapanFold ESMFold2-Fast structures. The rerun replaced those first 15 structures with JapanFold predictions using the same settings and seeds as the other 165. All 36 ordered inputs retained seeds 0–4 and both candidate perspectives for every structure. JapanFold ran single-sequence Fast mode with one sample per seed, ten requested recycles, and fifteen requested diffusion steps. The provider logs recorded eleven trunk passes and eleven diffusion steps in all 30 batches. [JapanFold ranking](results/esmfold2-japanfold-b2ar-six-ranking.csv) · [Score matrix CSV](results/esmfold2-japanfold-b2ar-six-score-matrix.csv) · [Score matrix PNG](results/esmfold2-japanfold-b2ar-six-score-matrix.png) · [Earlier mixed ranking](results/esmfold2-mixed-b2ar-six-ranking.csv)

| Nanobody | JapanFold ESMFold2-Fast rank | Earlier mixed rank | Archived AlphaFold2-multimer rank | Experimental rank |
|---|---:|---:|---:|---:|
| Nb6B9 | 1 | 1 | 1 | 1 |
| Nb6A9 | 2 | 2 | 3 | 2 |
| Nb6A5 | 3 | 3 | 2 | 3 |
| Nb6B11 | 4 | 5 | 5 | 4 |
| Nb80 | 5 | 4 | 6 | 6 |
| Nb6A10 | 6 | 6 | 4 | 5 |

Spearman ρ against experiment was 33/35 (0.943) for the JapanFold panel, 29/35 (0.829) for the mixed panel, and 31/35 (0.886) for archived AlphaFold2-multimer. The JapanFold panel reached the correlation cutoff and placed Nb6B9 first. Nb80 moved from fourth to fifth. The qualification gate requires sixth place. ESMFold2-Fast remains experimental. The ColabFold2 preview also retains its separate free-Colab validation record.

JapanFold charged **$0.015712** for the 15 reruns and **$0.194861** in total for 2,698.25 processor-seconds across 30 settled receipts. The estimate for a future run of 180 similar predictions is **about $0.20**. JapanFold bills the actual processor time to the JapanFold account. [Estimated costs and billing](nanobody-runner.md#estimated-esmfold2-usage-costs) · [JapanFold pricing](https://japanfold.aiand.com/docs/models-and-limits/)

## Publication and public Colab

APPRAISE `1.3.2` was published through the approved GitHub Trusted Publisher. The `v1.3.2` tag points to `06bc9e024417c9fbbebf6250dd3e1991336f072a`. The two public downloads matched the tested workflow files byte for byte. A fresh PyPI install reproduced both CPU examples and resume. [PyPI release](https://pypi.org/project/appraise/1.3.2/) · [GitHub release](https://github.com/xz-ding/APPRAISE/releases/tag/v1.3.2) · [Publishing run](https://github.com/xz-ding/APPRAISE/actions/runs/45892092792)

| Published file | SHA-256 |
|---|---|
| Wheel | `5c0f05c74aafbe2355c09f4c53bf5f28f86bc2e3fa3dffb3a02fa661ad97b40b` |
| Source archive | `178b92e967a4680fadfde15df48aa864512cc60b7caa6327421f5cf33f3b5786` |

The canonical README links to the guided public notebook merged at `fdb8cd14dab4393cceb8588fc8ff3013ec8c51cd`. In a fresh Colab CPU runtime on Python 3.13.18, **Run all** installed APPRAISE `1.3.2` from PyPI in 14.83 seconds. It scored ten published structures in 69 seconds and produced 20 measurements, a rank table, and a heatmap. Nb6B9 ranked first and Nb80 second. All three output hashes matched the run's result record. [Public notebook](https://colab.research.google.com/github/xz-ding/APPRAISE/blob/main/notebooks/Colab_APPRAISE.ipynb) · [Output hashes and runtime](implementation-state.json) · [Canonical CI](https://github.com/xz-ding/APPRAISE/actions/runs/45892804374)

Colab declined a free T4 allocation during the November 26 release check. The November 27 allocation completed the [peptide and miniprotein runtime checks](#peptide-and-miniprotein-runtime-checks).
