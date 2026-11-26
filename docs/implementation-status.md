# APPRAISE release candidate

The November 25 owner-supplied plan is recorded byte for byte in `AGENT_EXECUTION_PLAN.md`. Its SHA-256 is `3337349356d7503966539223616ed40017c56cf653a7e9ae355038ca650a66e2`. Work started from `8898c203d35b51d64b4766342ed9a9a32630cc13` on canonical `main`. The package and notebook are on the `codex/appraise-colab-runner` branch. [Draft PR](https://github.com/xz-ding/APPRAISE/pull/1) · [Machine-readable checkpoint](implementation-state.json)

## What runs today

The notebook, noninteractive command and Python interface share the published β₂AR workflow. Each sample records the producer, input role and copy, observed chain IDs, structure and measurement hashes. Recovery tests cover interruptions, damaged files and changed chain metadata. The preview adapter requests one ESMFold2-Fast diffusion sample per seed and the fallback requests one AlphaFold2 checkpoint per prediction. The guided test case contains three protein chains and 554 residues. [Candidate notebook](https://colab.research.google.com/github/xz-ding/APPRAISE/blob/0529084684d3f8a0c3fa08b8ae8d461703139146/notebooks/Colab_APPRAISE_candidate.ipynb) · [Usage](nanobody-runner.md)

| Evidence | Result |
|---|---|
| Published six-nanobody measurements | All 36 ordered inputs, 180 structures and 360 measurements represented. The APPRAISE 1.2 order is Nb6B9, Nb6A5, Nb6A9, Nb6A10, Nb6B11, Nb80. Spearman ρ against the experimental order is 31/35. |
| CPU structure-to-ranking run | Headless PyMOL measured the full chains from both directions in ten archived structures. All 20 perspectives were retained. Nb6B9 ranked first and Nb80 second. |
| Package | The wheel and source archive are built from a clean source export. The release check compares their package contents with source, installs the wheel outside the checkout, runs both examples and verifies resume. The hosted job retains the exact artifacts and their hashes. |
| Live free Colab | The built-in browser connected to a free Tesla T4 with 15 GiB VRAM and 12.67 GiB system RAM. Temporary checkpoint write and read-back passed. |
| AlphaFold2-multimer v2 | ColabFold 1.6.3 used templates and completed both input orders. It produced ten distinct 554-residue structures, 20 measurements, the ranked table and heatmap. Nb6B9 placed first. |
| ESMFold2-Fast | The v3.1.20 port with the streamed ESM-C loader completed both input orders. It produced ten distinct 554-residue structures, 20 measurements, the ranked table and heatmap. Nb6B9 placed first. |
| Six-nanobody comparison | JapanFold ESMFold2-Fast completed 36 ordered inputs, 180 distinct structures and 360 measurements. Spearman ρ was 33/35 against the experimental ranking. Nb80 finished fifth. The qualification gate requires sixth place. [Ranking](results/esmfold2-japanfold-b2ar-six-ranking.csv) |

All 40 scoped checks and five subtests passed locally and on GitHub’s Linux runner. They cover the two notebooks, backend I/O, the nanobody workflow, and archived measurement parsing. They also test ties, missing coverage, ensemble identities, checksum recovery, and signed colors in the guided and manual score matrices. [Hosted CPU checks](https://github.com/xz-ding/APPRAISE/actions/runs/45877547525)

## Free-Colab validation

The original port exhausted host memory while loading ESM-C. Streaming the same weight records resolved that stage. The original model graph expected 447 features from a 451-feature checkpoint. The pinned v3.1.20 port contains the upstream shape fix. The first new ESMFold2-Fast model and the first AlphaFold2-multimer model both passed written-chain and native-confidence checks. [Upstream checkpoint fix](https://github.com/sokrypton/alphafold3/commit/d87fa1eacaf275847f69e7f3ea854c362985de18)

The ESMFold2-Fast run took 58 minutes and 53 seconds after setup, and sampled peak GPU memory was 11,333 MiB. All ten structures passed chain and hash verification. The pinned environment was newly installed in the session and reused the model weights and ESM-C disk cache built earlier in that runtime. The table, measurements, heatmap and manifest hashes match the live notebook. [Ranking CSV](results/colab-t4-esmfold2-fast-b2ar-two-ranking.csv) · [Runtime and file hashes](implementation-state.json)

The AlphaFold2-multimer v2 test resumed from the first saved structure using ColabFold alignments and templates prepared for the exact two inputs. The resumed segment took 36 minutes and 18 seconds, excluding this preparation. Sampled peak GPU memory was 3,083 MiB. After Colab reconnected to the same runtime, the first structure's checksum still matched and all ten structures and final artifacts passed verification. [Ranking CSV](results/colab-t4-alphafold2-multimer-b2ar-two-ranking.csv) · [Runtime and file hashes](implementation-state.json)

The T4 handled the measured 554-residue three-chain workload. The guided ESMFold2 timing implied about 17 T4 GPU hours for the 180-prediction panel. Colab's free usage limit interrupted the run after 18 structures and 36 measurements.

The six-nanobody run used `/content/appraise-local/b2ar/b2ar-six-esmfold2-fast` on the GPU runtime. The automatic approval review declined the Drive filesystem mount because it requested account-wide file access. Six checksum-verified checkpoint archives were saved through the Drive web editor under `My Drive/APPRAISE/b2ar`. The fresh CPU runtime verified all 45 restored structure, confidence and measurement files across the 15 saved samples. Those samples form the durable T4 checkpoint.

On November 26, Colab again declined a free T4 allocation because of the account's GPU usage limit. The pinned APPRAISE `resume` command on a fresh Colab CPU runtime reused three ordered inputs, 15 structures and 30 measurements under the original run ID. All 45 file hashes and the 15 completed sample records still matched. The command exited at prediction 16 with `gpu_unavailable`.

For a T4-only replay, restore the six archives from `001-001` through `011-015` in ascending order into the same run directory. Each Markdown file starts with a marker and JSON metadata, followed by a base64-encoded gzip archive and an end marker. Verify each archive against the SHA-256 in the execution record before unpacking it. The last archive supplies the run manifest for all 15 saved samples. Install APPRAISE at `9c7e2e9590d90b02c46d2cf00e3be904867452ee` with the tested ESMFold2-Fast environment and run `python -m appraise resume --run-dir /content/appraise-local/b2ar/b2ar-six-esmfold2-fast`. The saved run ID is `1771fadf-bdf3-422a-bec7-c7c12b8980bd` and the next prediction is 16. [Archive checksums](implementation-state.json)

The guided runner covers the exact published nanobodies, their subsets and their verified receptor construct. Broader custom receptors and a measured peptide/miniprotein operating range remain open. Synthetic tests verify legacy ESMFold v1 export and recovery from exact source sequences. The manual notebook retains the existing advanced engines. The SKEMPI 2.0 follow-up and its within-target curation rules are in section 7 of the execution plan.

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

## Publication handoff

The candidate version is `1.3.2`. The workflow builds and tests the tagged revision, checks that it belongs to canonical main, verifies the artifact hashes and passes those files to PyPI using OIDC. The PyPI owner dashboard showed zero registered Trusted Publishers. The configuration awaiting confirmation uses owner `xz-ding`, repository `APPRAISE`, workflow filename `publish.yml`, environment `pypi`, and project `appraise`. [Publishing workflow](../.github/workflows/publish.yml) · [PyPI publisher settings](https://pypi.org/manage/project/appraise/settings/publishing/)

GitHub's `pypi` environment allows exactly the `v1.3.2` tag. PyPI and GitHub both returned HTTP 404 for the proposed version or tag before release. The release workflow also checks the canonical repository, version tag and merged source.

After the corrections release gates and publication access are verified, merge the package PR, tag the tested merged revision `v1.3.2` and dispatch `publish.yml` from that tag. Compare PyPI’s file digests with the workflow report and test a fresh install from PyPI. Activate the public notebook with the exact published version in a separate change and verify its canonical Colab link.
