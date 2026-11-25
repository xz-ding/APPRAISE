# APPRAISE release candidate

The November 25 owner-supplied plan is recorded byte for byte in `AGENT_EXECUTION_PLAN.md`. Its SHA-256 is `3337349356d7503966539223616ed40017c56cf653a7e9ae355038ca650a66e2`. Work started from `8898c203d35b51d64b4766342ed9a9a32630cc13` on canonical `main`. The package and notebook are on the `codex/appraise-colab-runner` branch. [Draft PR](https://github.com/xz-ding/APPRAISE/pull/1) · [Machine-readable checkpoint](implementation-state.json)

## What runs today

The notebook, noninteractive command and Python interface share the published β₂AR workflow. Each sample records the producer, input role and copy, observed chain IDs, structure and measurement hashes. Recovery tests cover interruptions, damaged files and changed chain metadata. The preview adapter requests one ESMFold2-Fast diffusion sample per seed and the fallback requests one AlphaFold2 checkpoint per prediction. The guided test case contains three protein chains and 554 residues. [Candidate notebook](https://colab.research.google.com/github/xz-ding/APPRAISE/blob/6a13198a8d09fea3573ab389d842c08b231cc3a9/notebooks/Colab_APPRAISE_candidate.ipynb) · [Usage](nanobody-runner.md)

| Evidence | Result |
|---|---|
| Published six-nanobody measurements | All 36 ordered inputs, 180 structures and 360 measurements represented. The APPRAISE 1.2 order is Nb6B9, Nb6A5, Nb6A9, Nb6A10, Nb6B11, Nb80. Spearman ρ against the experimental order is 31/35. |
| CPU structure-to-ranking run | Headless PyMOL measured the full chains from both directions in ten archived structures. All 20 perspectives were retained. Nb6B9 ranked first and Nb80 second. |
| Package | The wheel and source archive are built from a clean source export. The release check compares their package contents with source, installs the wheel outside the checkout, runs both examples and verifies resume. The hosted job retains the exact artifacts and their hashes. |
| Live free Colab | The candidate notebook renders in the built-in browser. Google sign-in and a free GPU allocation are pending. Drive persistence also awaits a live receipt. |
| New ESMFold2 scientific comparison | Zero newly folded benchmark structures. Hardware support, peak memory and ranking performance remain unmeasured. |

All 36 scoped checks and five subtests passed locally and on GitHub’s Linux runner. They cover the two notebooks, backend I/O, the nanobody workflow, and archived measurement parsing. They also test true ties, missing coverage, explicit ensemble identities and checksum recovery. The legacy capsid-ranking regressions were outside this validation. [Hosted CPU checks](https://github.com/xz-ding/APPRAISE/actions/runs/45721724297)

## Resume free-Colab validation

1. The account holder signs in to Google in the built-in browser. The agent then selects a zero-charge GPU runtime, verifies the allocation and performs the Drive write/read check.
2. Run the 554-residue pilot, verify native confidence, chain roles, memory and timing, then validate the two reciprocal jobs. Preserve the exact settings and all attempts under `MyDrive/APPRAISE/b2ar`.
3. Exercise the template-enabled AlphaFold2 fallback in a clean session. Freeze the ESMFold2 candidate before the 36-input scientific comparison. The acceptance threshold is 29/35 with Nb6B9 uniquely first and Nb80 uniquely last. Update hardware selection from the measured outcome.

The new guided runner covers the exact published nanobodies, their subsets and their verified receptor construct. Broader custom receptors, a measured peptide/miniprotein operating range, and legacy ESMFold-v1 boundary recovery remain outside the verified candidate. The manual notebook continues to expose existing advanced engines. The separate proposed SKEMPI 2.0 follow-up and its within-target curation rules are in section 7 of the execution plan.

## Publication handoff

The candidate version is `1.3.2`. The workflow builds and tests the tagged revision, checks that it belongs to canonical main, verifies the artifact hashes and passes those files to PyPI using OIDC. The PyPI Trusted Publisher fields are owner `xz-ding`, repository `APPRAISE`, workflow filename `publish.yml`, environment `pypi`, project `appraise`. The owner-side publisher configuration is unverified. [Publishing workflow](../.github/workflows/publish.yml) · [PyPI publisher settings](https://pypi.org/manage/project/appraise/settings/publishing/)

GitHub reported no deployment environments at this checkpoint. Automatic approval review declined environment provisioning because the creation request did not itself specify the release-tag pattern. The release workflow checks the canonical repository, version tag and merged source independently.

After qualification and publication access are verified, merge the package PR, tag the tested merged revision `v1.3.2` and dispatch `publish.yml` from that tag. Compare PyPI’s file digests with the workflow report and test a fresh install from PyPI. Activate the public notebook with the exact published version in a separate change and verify its canonical Colab link.
