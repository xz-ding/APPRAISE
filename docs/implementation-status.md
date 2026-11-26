# APPRAISE release candidate

The November 25 owner-supplied plan is recorded byte for byte in `AGENT_EXECUTION_PLAN.md`. Its SHA-256 is `3337349356d7503966539223616ed40017c56cf653a7e9ae355038ca650a66e2`. Work started from `8898c203d35b51d64b4766342ed9a9a32630cc13` on canonical `main`. The package and notebook are on the `codex/appraise-colab-runner` branch. [Draft PR](https://github.com/xz-ding/APPRAISE/pull/1) · [Machine-readable checkpoint](implementation-state.json)

## What runs today

The notebook, noninteractive command and Python interface share the published β₂AR workflow. Each sample records the producer, input role and copy, observed chain IDs, structure and measurement hashes. Recovery tests cover interruptions, damaged files and changed chain metadata. The preview adapter requests one ESMFold2-Fast diffusion sample per seed and the fallback requests one AlphaFold2 checkpoint per prediction. The guided test case contains three protein chains and 554 residues. [Candidate notebook](https://colab.research.google.com/github/xz-ding/APPRAISE/blob/c2c429fa2de7be3458d203a465da7e14436b5896/notebooks/Colab_APPRAISE_candidate.ipynb) · [Usage](nanobody-runner.md)

| Evidence | Result |
|---|---|
| Published six-nanobody measurements | All 36 ordered inputs, 180 structures and 360 measurements represented. The APPRAISE 1.2 order is Nb6B9, Nb6A5, Nb6A9, Nb6A10, Nb6B11, Nb80. Spearman ρ against the experimental order is 31/35. |
| CPU structure-to-ranking run | Headless PyMOL measured the full chains from both directions in ten archived structures. All 20 perspectives were retained. Nb6B9 ranked first and Nb80 second. |
| Package | The wheel and source archive are built from a clean source export. The release check compares their package contents with source, installs the wheel outside the checkout, runs both examples and verifies resume. The hosted job retains the exact artifacts and their hashes. |
| Live free Colab | The built-in browser connected to a free Tesla T4 with 15 GiB VRAM and 12.67 GiB system RAM. Temporary checkpoint write and read-back passed. |
| AlphaFold2-multimer v2 | ColabFold 1.6.3 used templates and completed both input orders. It produced ten distinct 554-residue structures, 20 measurements, the ranked table and heatmap. Nb6B9 placed first. |
| ESMFold2-Fast | The v3.1.20 port with the streamed ESM-C loader completed both input orders. It produced ten distinct 554-residue structures, 20 measurements, the ranked table and heatmap. Nb6B9 placed first. |
| New six-nanobody comparison | Colab interrupted the free T4 run after 18 structures and 36 measurements, then refused a GPU reconnect because of usage limits. The first 15 structures and 30 measurements were restored from Drive on a fresh CPU runtime. |

All 38 scoped checks and five subtests passed locally and on GitHub’s Linux runner. They cover the two notebooks, backend I/O, the nanobody workflow, and archived measurement parsing. They also test true ties, missing coverage, explicit ensemble identities and checksum recovery. The legacy capsid-ranking regressions were outside this validation. [Hosted CPU checks](https://github.com/xz-ding/APPRAISE/actions/runs/45728726304)

## Free-Colab validation

The original port exhausted host memory while loading ESM-C. Streaming the same weight records resolved that stage. The original model graph expected 447 features from a 451-feature checkpoint. The pinned v3.1.20 port contains the upstream shape fix. The first new ESMFold2-Fast model and the first AlphaFold2-multimer model both passed written-chain and native-confidence checks. [Upstream checkpoint fix](https://github.com/sokrypton/alphafold3/commit/d87fa1eacaf275847f69e7f3ea854c362985de18)

The ESMFold2-Fast run took 58 minutes and 53 seconds after setup, and sampled peak GPU memory was 11,333 MiB. All ten structures passed chain and hash verification. The pinned environment was newly installed in the session and reused the model weights and ESM-C disk cache built earlier in that runtime. The table, measurements, heatmap and manifest hashes match the live notebook. [Ranking CSV](results/colab-t4-esmfold2-fast-b2ar-two-ranking.csv) · [Runtime and file hashes](implementation-state.json)

The AlphaFold2-multimer v2 test resumed from the first saved structure using ColabFold alignments and templates prepared for the exact two inputs. The resumed segment took 36 minutes and 18 seconds, excluding this preparation. Sampled peak GPU memory was 3,083 MiB. After Colab reconnected to the same runtime, the first structure's checksum still matched and all ten structures and final artifacts passed verification. [Ranking CSV](results/colab-t4-alphafold2-multimer-b2ar-two-ranking.csv) · [Runtime and file hashes](implementation-state.json)

The 36-input ESMFold2 test is estimated at about 17 GPU hours from the guided test's median model time. It requires 180 complete predictions and a Spearman correlation of at least 29/35 with Nb6B9 uniquely first and Nb80 uniquely last. Default selection follows the measured configuration and scientific outcome.

The six-nanobody run used `/content/appraise-local/b2ar/b2ar-six-esmfold2-fast` on the GPU runtime. The automatic approval review declined the Drive filesystem mount because it requested account-wide file access. Six checksum-verified checkpoint archives were saved through the Drive web editor under `My Drive/APPRAISE/b2ar`. The fresh CPU runtime verified all 45 restored structure, confidence and measurement files across the 15 saved samples. The previous GPU runtime's working directory was absent. Predictions 16 through 18 will be recomputed.

On November 26, Colab again declined a free T4 allocation because of the account's GPU usage limit. The pinned APPRAISE `resume` command on a fresh Colab CPU runtime reused three ordered inputs, 15 structures and 30 measurements under the original run ID. All 45 file hashes and the 15 completed sample records still matched. The command exited at prediction 16 with `gpu_unavailable`.

When a free GPU is available, restore the six archives from `001-001` through `011-015` in ascending order into the same run directory. Each Markdown file starts with a marker and JSON metadata, followed by a base64-encoded gzip archive and an end marker. Verify each archive against the SHA-256 in the execution record before unpacking it. The last archive supplies the run manifest for all 15 saved samples. Install APPRAISE at `9c7e2e9590d90b02c46d2cf00e3be904867452ee` with the candidate notebook's pinned ESMFold2-Fast environment and run `python -m appraise resume --run-dir /content/appraise-local/b2ar/b2ar-six-esmfold2-fast`. The saved run ID is `1771fadf-bdf3-422a-bec7-c7c12b8980bd` and the next prediction is 16. [Archive checksums](implementation-state.json)

The new guided runner covers the exact published nanobodies, their subsets and their verified receptor construct. Broader custom receptors, a measured peptide/miniprotein operating range, and legacy ESMFold-v1 boundary recovery remain outside the verified candidate. The manual notebook continues to expose existing advanced engines. The separate proposed SKEMPI 2.0 follow-up and its within-target curation rules are in section 7 of the execution plan.

## Publication handoff

The candidate version is `1.3.2`. The workflow builds and tests the tagged revision, checks that it belongs to canonical main, verifies the artifact hashes and passes those files to PyPI using OIDC. The PyPI Trusted Publisher fields are owner `xz-ding`, repository `APPRAISE`, workflow filename `publish.yml`, environment `pypi`, project `appraise`. The owner-side publisher configuration is unverified. [Publishing workflow](../.github/workflows/publish.yml) · [PyPI publisher settings](https://pypi.org/manage/project/appraise/settings/publishing/)

GitHub reported no deployment environments at this checkpoint. Automatic approval review declined environment provisioning because the creation request did not itself specify the release-tag pattern. The release workflow checks the canonical repository, version tag and merged source independently.

After qualification and publication access are verified, merge the package PR, tag the tested merged revision `v1.3.2` and dispatch `publish.yml` from that tag. Compare PyPI’s file digests with the workflow report and test a fresh install from PyPI. Activate the public notebook with the exact published version in a separate change and verify its canonical Colab link.
