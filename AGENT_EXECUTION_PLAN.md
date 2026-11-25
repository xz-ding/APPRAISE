# APPRAISE and Colab-APPRAISE execution plan

Revised November 25, 2026 following the owner's clarification decisions.

**Canonical repository:** [xz-ding/APPRAISE](https://github.com/xz-ding/APPRAISE), default branch `main`. The reviewed head is `8898c203d35b51d64b4766342ed9a9a32630cc13`; discover and record the actual starting head when execution begins.

**Document status:** Specification for the implementation agent. ESMFold2 compatibility, benchmark results, operating limits, and release checks remain implementation deliverables. Values describing the historical AlphaFold2 reference below are reference expectations to verify, not new ESMFold2 measurements.

## 0. Agreed outcomes and scope

1. **Repair the public workflow.** A newcomer can open Colab-APPRAISE, complete setup, run a small supplied example, and find a readable ranked table and heatmap. The normal path runs one modeling engine. Advanced engines stay inactive unless selected, including under Run all.
2. **Use the six-nanobody panel as the key qualification.** Retain the published activated β₂-adrenergic receptor panel and the initial acceptance criterion in section 2. It is sufficient for this release's initial recommendation decision. Describe its biological scope explicitly.
3. **Optimize for small protein competition panels.** The initial workflow covers peptides, miniproteins, and nanobodies against defined receptor constructs. Establish and document the operating range by measurement. Large antibodies, multichain-per-candidate assemblies, large pooled screens, and arbitrary receptor sizes are not new qualification requirements.
4. **Evaluate ESMFold2-Fast as the preferred default candidate.** Promote one maintained integration only after the guided workflow and scientific gate pass, and only on its qualified runtime/input configuration. Retain AlphaFold2-multimer through ColabFold as the fallback where its own workflow is supported. CPU sessions can quantify and score saved structures.
5. **Release from the personal repository.** Code, PRs, release tags, public Colab links, and the package publication workflow originate from `xz-ding/APPRAISE`. The laboratory fork is an optional downstream copy; synchronizing it is not a release prerequisite.
6. **Support agent users directly.** Ship one small runner shared by Colab, Python, and a simple command, with explicit inputs/settings, resumable execution, structured results, and clear failure status. Defer a general CLI framework and a new YAML configuration system.
7. **Proceed autonomously after initial access setup.** The implementation agent owns implementation, checks, PRs, merge, publication, and public verification. Routine decisions require no additional sign-off. Account interactions and repository-enforced approvals have narrowly defined handoffs.
8. **Keep GPU spending at $0.** Use genuine free Google Colab for the required GPU evidence. No paid fallback is authorized. Complete independent work during external interruptions, preserve checkpoints, and use the finite blocked-state policy below.
9. **Recommend a larger benchmark without making it a release dependency.** Use a curated subset of SKEMPI 2.0 for a subsequent broader affinity evaluation. Section 7 specifies the proposal and its limitations.

### Meaning of a recommendation

Distinguish operational support from scientific validation. Peptide or miniprotein inputs may pass runtime and format checks without having been scientifically benchmarked for affinity ranking.

Use concise wording such as:

> Recommended engine on this supported configuration: ESMFold2-Fast. Initial ranking qualification: six related nanobodies against activated β₂AR. Performance on other targets and binder classes remains to be characterized.

APPRAISE reports relative binding propensity within a competition panel. It does not produce a measured or calibrated absolute affinity. A top rank alone does not establish binding without appropriate controls. See the [APPRAISE paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11184338/).

## 1. Establish access and execution ownership first

### 1.1 Access inventory

At execution start, inspect the actual checkout, preserve existing user edits, and establish the following usable routes before long GPU work or publication.

| Dependency | Establish early | Evidence to retain |
|---|---|---|
| Canonical GitHub repository | Access to `xz-ding/APPRAISE`, feature branch and PR creation, actual merge requirements, and available CI. | Starting commit, branch, PR URL when created, and required checks/approvals. |
| PyPI project `appraise` | A working owner-authorized publication route. Prefer the project's configured Trusted Publisher; prepare a suitable workflow if setup is needed. | Workflow/repository/environment identity and whether owner-only setup remains. Record no secrets. |
| Free Colab | An authenticated notebook session through the signed-in product or an officially supported connection, and the ability to request a free GPU. | Actual connection route, assigned runtime, and any precise account action still needed. A local notebook is not free-Colab evidence. |
| Persistent checkpoints | A specific writable destination, normally a dedicated Drive project folder, accessible from the notebook after runtime loss. | A small write/read check, the durable path, and how to locate and resume the run. |

The connected GitHub permissions observed during planning allow writing to the personal repository. Recheck at execution time. PyPI and Colab access are not established by this document.

Prepare all owner-reviewable configuration before requesting an account action. Consolidate missing setup into one short handoff whenever possible. Never ask for passwords or publishing tokens in chat. Keep credentials, downloaded weights, and large model outputs outside version control.

Continue independent implementation, CPU checks, documentation, and packaging while an access dependency is unresolved. Access discovery must happen early; an external dependency must not unnecessarily block independent work.

### 1.2 Durable execution record

Maintain a short human-readable progress record and a machine-readable execution state in persistent storage. Include the current tested commit, PR, completed phases, CPU results, Colab evidence, benchmark status, selected defaults, checkpoint paths, blockers, and exact next action. Small non-sensitive reports can be committed or linked from the PR; large structures belong in durable artifact storage.

Distinguish the project execution record from the per-run manifests in section 4. Both must survive loss of a Colab VM or an agent context.

### 1.3 Interruption and blocker policy

| State | Meaning | Required action |
|---|---|---|
| Running | Prerequisites are available and useful work can proceed. | Continue without discretionary approval pauses. |
| Recoverably interrupted | A transient failure occurred; saved work and an immediately usable recovery route exist. | Diagnose, resume from validated checkpoints, and retain job identities. |
| Blocked by an external prerequisite | Authentication, unavailable allocation, publishing access, or an enforced approval prevents the next dependent step. | Name the dependency, action owner, checkpoint, and exact resume action; continue independent work. |
| ESMFold2 evaluated but ineligible | A reproducible runtime limitation or the declared scientific gate failed. | Report the measured result and retain a qualified fallback. |
| Complete | All outcomes declared for the particular release were verified publicly. | Report public versions, links, evidence, and any separately unresolved qualification. |

Use a bounded transient retry policy: one diagnosed retry after the initial failure for an unchanged access/allocation attempt, or at most two retries for a modeling job when there is a concrete recoverable cause. Record attempts. A changed external condition can justify a new attempt; an unchanged quota denial cannot justify endless polling.

When all remaining work depends on an unavailable service or owner action, stop active retries and issue one precise handoff. The checkpoint must identify what is complete and how to resume. Do not promise automatic future resumption unless a durable execution mechanism has actually been configured and tested.

There is no overall project deadline or scheduled pause. This does not imply uninterrupted Colab allocation. Follow [Colab's resource and session policies](https://research.google.com/colaboratory/faq.html); preserve work rather than attempting to bypass them.

A corrections release may finish with ESMFold2 qualification explicitly blocked or pending, provided its own public default workflow passes. Publication blocked by missing credentials is reported as publication blocked, even if implementation and validation are complete.

## 2. Reproduce the scientific reference and freeze qualification

### 2.1 Historical reference

Use the **six nanobodies binding activated β₂AR**, Figure 3H–J, throughout the initial qualification. A two-candidate subset supplies the guided example; the complete panel determines whether to recommend the new engine.

| Rank, strongest first | Experimental order | Published AlphaFold2-multimer + APPRAISE 1.2 order |
|---:|---|---|
| 1 | Nb6B9 | Nb6B9 |
| 2 | Nb6A9 | Nb6A5 |
| 3 | Nb6A5 | Nb6A9 |
| 4 | Nb6B11 | Nb6A10 |
| 5 | Nb6A10 | Nb6B11 |
| 6 | Nb80 | Nb80 |

The paper reports Spearman's ρ = 0.89; these orders give exactly **31/35 = 0.885714…**. Experimental ranks are the bracketed labels in Figure 3H. The manuscript input files provide identities and sequences. [Paper and Figure 3](https://pmc.ncbi.nlm.nih.gov/articles/PMC11184338/#fig3) · [Published figure](https://cdn.ncbi.nlm.nih.gov/pmc/blobs/f91d/11184338/15dcfc94dae3/gr3.jpg)

The supplied plan identifies a reusable measurement CSV from Figure 3 archive part 2 with:

- 180 distinct model identities;
- two candidate perspectives per structure, totaling 360 measurements;
- 36 ordered candidate–competitor groups, with 10 measurement rows per group after the perspective convention;
- SHA-256 `9c5286291c2a36a0a5096ad7090baf42f5b4fb2fb857e631cb1e116d1ce7221e`.

Fetch the [individual reference CSV](https://data.caltech.edu/api/records/kxjgj-tfk18/files/Figure_3_related%28non-AAV_protein_applications%29_part_2_of_3.zip/container/b2AR_6x6/database_APPRAISE_measurements_09271324_note_N-terminal_anchor_site.csv), verify its checksum and inventory, and reproduce the published order through the package's normal scoring entry point. Save a CPU regression test. [CaltechDATA record](https://data.caltech.edu/records/kxjgj-tfk18)

The archived AlphaFold2 panel is the scientific comparator; refolding that complete panel is unnecessary for the initial release. Fresh AlphaFold2 folding is still required to validate the repaired public fallback workflow.

### 2.2 Accepted initial qualification criterion

Freeze the following before inspecting ESMFold2 rankings:

- All six candidates have the required complete input/model coverage.
- Nb6B9 is uniquely first and Nb80 uniquely last.
- Spearman's **ρ ≥ 29/35 ≈ 0.828571** against the experimental order.
- Use unrounded values for the decision and average ranks for genuine unresolved ties.
- Report the experimental, historical AlphaFold2, and ESMFold2 orders; both correlations; their difference; ties; and structural/runtime failures.

This is the agreed initial qualification rule. It permits one attainable increment below the published untied six-item correlation. It does not establish statistical equivalence to AlphaFold2 or generalization across targets. Six related binders against one receptor remain one biological panel, regardless of the number of modeled structures.

A scientific miss must remain a reported result. Do not relax the threshold, change the panel, or tune settings against the known ordering after viewing results. Broader affinity benchmarking is a separate stage.

### 2.3 Exact benchmark inputs and sampling

Use the [36 manuscript FASTAs](https://github.com/xz-ding/APPRAISE/tree/main/demo/example_input_sequences_from_manuscript/b2AR_6x6/b2AR_6x6_input_fasta) and preserve their chain order. Each complex contains two 120-residue nanobodies and the 314-residue receptor construct: **554 residues in three chains**.

Use the manuscript receptor settings: N-terminal anchor, Dmax 93.8 Å, axial ratio 2.23, and Rminor 21 Å. Match the exact construct sequence, not the receptor name alone. [Receptor settings](https://github.com/xz-ding/APPRAISE/blob/main/demo/example_input_sequences_from_manuscript/receptor_sequences_part_1.csv)

Keep APPRAISE 1.2 equations, reciprocal inputs, and the paper's self-pair convention. The AlphaFold2 archive contains five models per ordered input.

For ESMFold2, predeclare five distinct seeds and one diffusion sample per seed per ordered input. If the chosen upstream API uses a different sampling convention, establish an equivalent five-output protocol before freezing it. Verify that wrappers retain all requested samples rather than silently selecting a best-of-many result. Execute sequentially when required for memory.

Record the exact model variant, implementation, structure and language-model checkpoint revisions, precision, diffusion steps, refinement loops, and all consequential settings. A model identity consists of its ordered input, protocol fingerprint, seed, and sample index. Retries retain this identity.

Require **36 × 5 = 180 structures and 360 candidate-perspective measurements** before applying the scientific gate. Mirrored candidate rows share one structure identity; they are not independent model replicates. Record reciprocal comparisons and self-pairs explicitly when computing counts.

A smoke-test prediction can be reused only if its input and complete frozen protocol match. Do not retain favorable pilot outputs while excluding other valid outputs from the declared protocol.

### 2.4 Correctness and historical compatibility

Reproduce the reference before changing scoring-related behavior. Use synthetic, independently interpretable examples to justify correctness fixes.

Investigate an unexpected change in the historical order before updating expected results. If a proven correction necessarily changes historical behavior, preserve a labeled reproduction of the original behavior and document the corrected result and cause. Do not tune parameters to recover the expected order. This does not authorize a new engine-specific scoring model or relaxation of the qualification rule.

## 3. Repair the existing package and public notebook

Start on a feature branch from the current canonical target branch. Record the baseline test outcome and open a draft PR early. Preserve existing user edits and demo data. Use scoped commits and record their IDs in progress updates or the PR.

| Work | Required change | Evidence required |
|---|---|---|
| Notebook/package bootstrap | Candidate runs install the exact candidate wheel or source revision and show installed version/source. Package the guided inputs, saved example, and verified receptor settings so a wheel-only installation works outside the checkout. Public notebook installs its corresponding exact published package version. | Clean environment outside the checkout imports a candidate-specific feature and loads the examples/settings. Repeat using the published package. |
| ColabFold installation | Pin a tested release or immutable revision. Resolve template tools and optional Amber dependencies for actual Colab Python/CUDA. Keep the default installer from replacing Colab's Python. Preserve legacy model choices under advanced settings. | A clean free-Colab templates-enabled AlphaFold2-multimer run reaches a saved structure and ranking. Exercise changed optional install branches and label their actual validation status. |
| Notebook execution | Replace position-sensitive cell magics, define the input-reader logger before use, and preserve seeds/inputs on repeat Step 1 execution. Guard every modeling block by the selected engine. Locate cells by content or stable metadata and avoid unrelated notebook reserialization. | Run changed cell bodies/shared functions, check defaults and inactive branches, verify notebook round trips, and include an input directory containing an unrelated file. |
| Chain identity | Carry ordered chain IDs, names, roles, exact sequences, and copies from job generation to scoring. Map output IDs to input roles even for reordered chains and identical self-pair sequences. | PDB/mmCIF fixtures cover reordering, duplicate sequence identities, missing metadata, and ambiguous mappings. Never infer a duplicate chain's identity from sequence alone. |
| Legacy ESMFold boundaries | For joined v1 outputs, remove only known artificial linker residues on export. Recover old boundaries from original separated inputs; mark ambiguous legacy cases. Retain published v1 inference settings. | Fixtures include terminal glycines, internal glycine runs, artificial linkers, and absent boundary metadata. |
| Confidence | Normalize once using the engine and writer's documented units and preserve provenance. Do not interpret genuine experimental B-factors as confidence. | Equivalent 0–1 and 0–100 fixtures give equivalent scores; experimental B-factors remain unchanged. |
| Scores and CSV schema | Resolve duplicated angle headers by matching values to emitted positions; retain old CSV readability and every competitor in pooled inputs. Count actual structures/perspectives before aggregation. Correct the tie-break interval. | Synthetic pooled inputs, incomplete grids, nonzero tie-point cases, and the archived nanobody reference establish behavior. |
| Outputs and failures | Fix one-column plotting, empty-input failures, broken example links, mismatched receptor diagnostics, and swallowed subprocess failures. | Small one/multiple-column plots and representative failure cases show the failing step, relevant message, and next action. |
| Shared runner | Move orchestration into package functions called by Colab and the thin Python/command interfaces in section 4. | CPU orchestration fixtures and a real Colab run exercise the same entry point. |

Use small synthetic structures for format/chain regressions and the archived β₂AR measurements for numerical regression. Run the existing tests, targeted new checks, and headless PyMOL quantification through the shared entry point in CI. Record which checks ran on CPU, in a local environment, and in actual Colab.

Preserve the advanced paths for AlphaFold2-multimer, legacy ESMFold, Boltz-1, Chai-1, and the existing external OpenFold3 handoff. Retain their input/output compatibility and label their support levels accurately. Do not imply that preserving an adapter qualifies its upstream model on current hardware.

## 4. Ship a small runner shared with Colab

### 4.1 One execution path

Implement one package pipeline for validation, input generation, modeling, output validation, quantification, ranking, and result reporting. Notebook cells collect inputs and present status/results through that pipeline.

Provide one typed settings object, one Python entry point, a thin command entry point, and a saved-structure scoring mode. The following is the proposed public interface to implement and document:

```python
from appraise.runner import RunConfig, run

config = RunConfig(
    example="b2ar-two",
    engine="auto",
    output_dir="runs/b2ar-demo",
    predictions_per_input=5,
    seed=0,
)
result = run(config)
```

```bash
python -m appraise run --example b2ar-two --engine auto --predictions-per-input 5 --output runs/b2ar-demo
python -m appraise resume --run-dir runs/b2ar-demo
python -m appraise score --example b2ar-saved --output runs/b2ar-scored
```

Finalize exact names during implementation and keep documentation and tests consistent. A console alias may call the same implementation. No separate service, plugin framework, general scheduler, or YAML system is required.

Define `predictions_per_input` as the total retained predictions per ordered input, not a multiplier applied separately to model checkpoints and seeds. For the five-prediction ESMFold2 protocol, base seed b expands to the declared seeds b through b+4, with one diffusion sample each. The AlphaFold2 adapter must resolve its model/seed ensemble into an explicit output list and total count; do not accidentally request five models for each of five seeds. Record model, seed, and sample identity for every output.

The relative paths above illustrate local invocation. In Colab, resolve the run directory under the verified persistent destination, or synchronize each validated job and updated manifest there immediately. If modeling uses local scratch, record both paths and verify the durable copy before marking checkpoint persistence complete. Preflight must show the destination that will actually receive the current run's checkpoints.

For users' own inputs, support existing candidate/receptor formats, exact receptor construct, output directory, selected engine, explicit required geometry/anchor parameters or verified preset, relevant model settings, and seed/sample settings. Provide one complete custom-input example. Reject missing scientific inputs with an actionable message rather than asking interactive command-line questions.

Resolve defaults into one effective configuration. Record all effective values; shared values must not diverge among notebook cells, API calls, and command invocations. `auto` resolves to an eligible engine once per run, and the selected engine/reason is displayed before modeling.

The minimal new runner must cover the guided default, qualified AlphaFold2 fallback, and saved-structure scoring. Existing advanced engine workflows may retain their documented manual paths for this release; do not expose an unimplemented command option merely to match their notebook labels.

### 4.2 Stable output contract

Return a structured result with at least:

- schema version and run ID;
- status and selected engine;
- expected/completed input jobs and model predictions;
- paths to measurements, rankings, heatmap, manifest, and logs;
- failure code and next action when applicable.

Write a run manifest with normalized inputs and hashes, ordered chain roles/copies, effective scientific settings, source/package version, engine implementation/checkpoints, sample identities, environment for each execution segment, and per-job completion states.

Recommended files are `run_manifest.json`, `result.json`, `measurements.csv`, `ranking.csv`, a heatmap, and a readable log. Preserve established CSV compatibility or document necessary additive schema changes. Include schema versions for machine-readable files.

Use the same result information to show human-readable status and links in Colab. The user should see which stage is running, where outputs are saved, and whether coverage is complete.

### 4.3 Resume behavior

A completed model can be reused only when its input/protocol fingerprint matches and its expected outputs pass integrity, chain/sequence, and metadata checks. File existence alone is insufficient.

Include relevant stage-specific implementation and format revisions in completion checks. An inference or serialization correction requires regenerating affected predictions or exports; a quantification/scoring correction requires invalidating and recomputing affected derived results from still-valid structures. Preserve unaffected folding work. Recording a new package version without invalidating stale derived results is insufficient.

Write files and completion records atomically where possible. Mark a model complete only after validation. On restart, reuse complete valid jobs and retry incomplete/interrupted jobs under their original identities. A changed engine, checkpoint, input, or scientific protocol requires a separate run. Record changed runtime environments even when the scientific protocol remains compatible.

Quantification must deduplicate or upsert by structure identity plus candidate perspective. Re-running a cell, retrying a model, or resuming a session cannot add extra evidence.

Persist checkpoints incrementally. Do not wait until the entire folding stage ends before saving results. Never present an incomplete panel as a completed qualification. Preserve partial measurements and report missing coverage; the six-nanobody acceptance gate requires full coverage.

### 4.4 Failures and focused tests

Use a small status set: `complete`, `interrupted`, `blocked`, and `failed`. The command exits zero only for the documented successful outcome. Document nonzero failure categories and include the precise reason in the structured result.

A CPU fixture or fake modeling adapter should verify:

1. Python and command interfaces exercise the shared pipeline and return equivalent results.
2. An interrupted run resumes without recomputing validated jobs or duplicating measurements.
3. Changed inputs/settings cannot reuse incompatible outputs.
4. Missing or malformed model files are not accepted as completed predictions.
5. Invalid inputs, missing dependencies, and unavailable resources produce useful non-success results.
6. Saved-structure quantification/scoring succeeds without a GPU.

These tests verify orchestration. Live Colab validation remains necessary for actual modeling and installation.

## 5. Qualify ESMFold2 and document the operating range

### 5.1 Select one supported implementation

Evaluate an immutable revision of the [ColabFold2 preview](https://github.com/sokrypton/ColabFold/blob/main/ColabFold2_preview.ipynb) ESMFold2-Fast integration first. If it cannot meet input, memory, confidence, or output requirements, evaluate [Biohub's documented implementation](https://github.com/Biohub/esm) and [local GPU tutorial](https://github.com/Biohub/esm/blob/main/cookbook/tutorials/esmfold2_local_gpu.ipynb).

ColabFold2's route is a port. Verify its actual implementation and weight sources, ESMC features, Fast checkpoint, confidence head, native multichain inputs, and output serialization. A selectable model name is not evidence that these components work together.

Use supported model loading and serialization; do not reimplement the neural model. Maintain one ESMFold2 route for this release. Label the precise variant and route in the notebook, manifest, and results.

Full ESMFold2 and ESMFold2-Fast are separate configurations. Fast is a single-sequence model; do not imply that it uses MSA conditioning, or silently replace it with a smaller variant missing required confidence outputs. [Official model card](https://huggingface.co/biohub/ESMFold2-Fast)

### 5.2 Measure on genuine free Colab

Use the GPU actually assigned. Record GPU model, usable VRAM, system RAM, Python/CUDA, exact dependency revisions, checkpoints, installation/download time, peak GPU memory, and inference time.

Start with one 554-residue complex and one declared sample. Try documented compatible low-memory settings if needed, then freeze the candidate protocol before the scientific comparison. Record failed inputs and stages as well as successful settings. Treat the full-model Hugging Face memory guidance as context, not a measured limit for Fast. [Hugging Face example](https://github.com/huggingface/notebooks/blob/main/examples/protein_folding.ipynb)

Do not claim that a local notebook, imported Colab package, or non-Colab GPU proves the public free-Colab path.

### 5.3 Prove the guided workflow

From a fresh Colab session, open the candidate notebook, install the exact candidate package, select the supplied example, generate competition inputs, fold them, quantify them, and save a linked ranked table and heatmap.

Use Nb6B9 and Nb80 with two reciprocal non-self inputs and five predictions per input: **10 distinct structures and 20 candidate-perspective measurements**. After reciprocal aggregation, each directed candidate–competitor group contains observations from 10 distinct structures, not five or twenty. Treat self-pairs as deliberately unscheduled, show their heatmap cells as unmeasured, and exclude them from completeness and match-percentage denominators. Verify counts and ranking behavior with a CPU fixture. If the validated shared scorer requires self-pairs, use the complete 2 × 2 panel and document **20 structures and 40 measurements**. The six-candidate benchmark always retains all 36 inputs and its 36 × 5 protocol.

Check written sequences and chain roles, confidence units, engine label, output paths, and completeness. Save the executed notebook or log with output hashes. A two-candidate example demonstrates operation and interpretation; it does not replace the six-candidate scientific gate.

Also provide a small saved-structure example that can quantify and score on CPU after installation.

### 5.4 Run the frozen six-nanobody benchmark

Run the section 2 protocol once for the selected candidate configuration. Checkpoint completed outputs to the persistent destination and resume across free sessions if necessary.

Classify problems as installation/orchestration, structure serialization/identity, reproducible runtime limits, or scientific ranking results. Keep all measured outcomes. A verified implementation defect can justify a documented correction and rerun; changing scientific parameters after seeing rankings creates a new, explicitly exploratory analysis.

### 5.5 Establish a conservative operating range

Alongside the 554-residue benchmark, use a few representative protein-only inputs covering the intended peptide/miniprotein/nanobody workflow and its proposed size range. These are runtime/format checks, not new affinity benchmarks. Do not perform an exhaustive size sweep.

Document the tested configurations in a concise support table:

| Field | Required record |
|---|---|
| Runtime | Actual free-Colab GPU, usable memory, Python/CUDA, dependency revisions. |
| Model | Engine variant, implementation, checkpoint, precision, and frozen inference settings. |
| Input | Exact sequences or hashes, receptor/binder lengths, total residues or tokens, and chain count. |
| Workload | Number of candidate pairs, predictions, and samples per input. |
| Resources | Installation time, per-prediction time, peak GPU memory, and relevant system-memory limit. |
| Outcome | Completed, recoverable failure, reproducible limit, or untested. |
| Interpretation | Operationally supported range and the narrower scope of scientific qualification. |

A largest successful input does not establish a universal maximum. Hardware memory alone does not establish compatibility. Candidate count mainly affects total runtime and persistence; individual complex size affects whether a prediction fits.

Before launching a screen, show exact input and prediction counts. For an ordered square panel with self-pairs, N candidates and S samples per input require **N² × S predictions**. Estimate completion time only from observed, relevant timings and label it as an estimate.

### 5.6 Choose defaults from evidence

Recommend the qualified Fast configuration automatically only within its recorded runtime/input range. Display the engine and selection reason before downloads or inference. A resumed run retains its engine and protocol.

On other configurations, use AlphaFold2-multimer only where its own workflow is supported. Do not assume it will fit an input merely because Fast failed. For untested/oversized inputs, explain the measured limitation and offer supported alternatives. Never silently truncate sequences, reduce sample count, or change precision to make an input pass.

A fallback starts a separate run and preserves the first run's evidence. Do not combine engines within a ranking. With no GPU, offer saved-structure scoring and instructions for starting a GPU session.

ESMFold2 outcomes are:

- **Qualified:** Guided workflow and scientific gate passed on the recorded configuration.
- **Evaluated but ineligible:** A reproducible runtime limit or scientific miss occurred; retain the qualified fallback and describe the result.
- **Blocked before qualification:** Access or allocation prevented a required test; compatibility remains unresolved.

A functioning but unqualified integration may ship as an advanced experimental option with its measured status. Higher-memory tiers remain unqualified until actually exercised at zero charge.

## 6. Ranked ten UX improvements

Rank by the effect on a first-time user or agent reaching a correct, understandable result. Ship the focused versions of **1–9** and the exact job-count component of **10**. Original item numbers preserve traceability to the earlier plan.

| Rank | Original item | Release scope | User benefit |
|---:|---|---|---|
| 1 | 3.9 · Guided example | Ship | One clear action after setup runs a small competition and ends with linked rankings and a heatmap. Inactive engines remain inactive under Run all. |
| 2 | 3.1 · Preflight | Ship | Before long downloads, show hardware, selected engine, measured support range, persistent destination, and the first corrective action. Offer CPU saved-structure scoring. |
| 3 | 3.10 · Quickstart and custom run | Ship | Put the canonical Open in Colab link, expected example result, and output location first. Provide a short complete path for the user's own inputs. |
| 4 | 3.8 · Shared settings | Ship | Inputs, engine, output directory, and common settings resolve consistently in Colab, Python, and the command. Model-specific settings have distinct names. |
| 5 | 3.3 · Durable resume | Ship | Reuse validated work after disconnects, recover partial jobs, and prevent duplicate model/measurement counts. |
| 6 | 3.6 · Small agent runner | Ship; promoted from deferred | One documented noninteractive command uses the shared pipeline, accepts explicit settings, resumes runs, and returns structured status/results. Defer a general CLI/YAML framework. |
| 7 | 3.5 · Actionable errors | Ship | Identify the failed stage, affected input, preserved work, and next step, including memory errors and insufficient coverage. |
| 8 | 3.7 · Receptor settings | Ship lookup and guide | Fill verified geometry by exact construct identity, show units, and explain how to supply defensible custom settings. Defer a general geometry estimator. |
| 9 | 3.4 · Progress and interpretation | Ship | Show active stage, completed/total predictions, elapsed time, and output links. Explain ties, missing comparisons, and the relative meaning of ranks. |
| 10 | 3.2 · Workload estimate | Ship exact counts; defer broad estimator | Show jobs/predictions before execution. Add timing estimates only when measurements support them. Defer general cost and Drive-capacity forecasting. |

For a custom receptor, explain each geometry/anchor parameter and how it maps to the modeled construct. Do not silently use a known receptor's settings for an unknown construct. Provide a worked example rather than requiring users to infer the settings from source code.

Keep advanced manual engines accessible without adding steps to the guided path. Defer archive-based batch storage, wholesale PyMOL performance work, additional confidence-based ranking models, and general scientific score redesign.

## 7. Larger affinity dataset recommendation: SKEMPI 2.0

### 7.1 Recommendation and scope

Use **a curated set of within-target panels from SKEMPI 2.0** for the next broader validation stage. It is an established benchmark source for mutation-induced changes in protein–protein affinity.

The original paper reports **7,085 records, 6,187 unique mutational entries, 237 interactions, and 345 PDB entries**. These are not counts of independent binders or ready-made APPRAISE panels. The database provides mutant/reference affinities and experimental metadata. [Original paper](https://doi.org/10.1093/bioinformatics/bty635) · [Official database](https://life.bsc.es/pid/skempi2/default/index) · [Field definitions](https://life.bsc.es/pid/skempi2/info/faq_and_help)

**AB-Bind** is an antibody-focused alternative: **1,101 mutants across 32 complexes** in its original publication. It overlaps SKEMPI and must not be counted automatically as an independent validation dataset. [AB-Bind paper](https://pubmed.ncbi.nlm.nih.gov/26473627/) · [Author-maintained data](https://github.com/sarahsirin/AB-Bind-Database)

The first release includes this recommendation and a clear follow-up protocol. It does not require folding the full database, preparing every eligible panel, or running a new large benchmark.

### 7.2 Construct valid APPRAISE panels

The following are proposed evaluation rules to freeze before a later modeling campaign:

1. **Group by exact biological and experimental context.** Use an unchanged target construct, a common binder parent/background and compatible binding site, and comparable assay/condition blocks. Target mutations create a separate target panel. Prefer binder-side variants of one parent where physical competition is plausible.
2. **Reconstruct exact sequences.** Resolve PDB and cleaned residue numbering, all chain roles, and each mutation. Verify the actual reference sequence; the database reference can itself be engineered. Retain source row IDs and structure provenance.
3. **Use affinity labels correctly.** Prefer pKD = −log10(KD in molar units), where higher is stronger. Within a shared reference and experimental context, −ΔΔG can define the same direction. ΔΔG = RT ln(KD_mutant / KD_reference). Do not compare ΔΔG across unrelated parents as absolute affinity, and do not label APPRAISE scores as predicted KD.
4. **Handle uncertainty and duplicates.** Prespecify treatment of duplicate sequences/conditions, experimental replicates, alternative PDB representations, and forward/reverse mutation records. Retain inequality/no-binding labels. Do not assign arbitrary exact affinities to censored values.
5. **Respect the competition assumption.** Do not combine unrelated binders against different epitopes into an affinity tournament merely because they share a target name. Confirm that candidate variants plausibly compete for the same site.
6. **Use biologically appropriate geometry.** Many targets are soluble. Use only a justified, supported scoring configuration; do not invent a membrane anchor. If suitable scoring would require a new scientific method, exclude that panel from the initial subset and document the reason.
7. **Separate model input from evaluation labels.** Freeze sequences, panel selection, sampling, scoring, and model/template settings before inspecting predictions. Use source structures for sequence/context verification without unintentionally granting one engine a privileged template.
8. **Describe historical-data limitations.** Record known training/template cutoffs and possible overlap. Treat these panels as historical validation, not proof of prospective generalization. Do not tune on a family and count close relatives as independent held-out evidence.

Freeze the source version/download hash, curation rules, and exclusions. Report the actual usable panel/variant counts after filtering. No number of eligible panels is promised before this audit.

### 7.3 Evaluate at the panel level

Compare qualified ESMFold2 and AlphaFold2 workflows on matched panels and a declared compute/sampling protocol. The historical β₂AR CSV covers only β₂AR; new panels need their own comparator predictions or independently suitable archived predictions.

Report per-panel Spearman correlation on exact-affinity labels, tie-aware pairwise agreement, top-binder recovery, completeness/failure rates, peak memory, and wall time. Censored labels can support comparisons only where their order is unambiguous. Prespecify how ties and uncertain experimental ordering are handled.

Summarize across panels with equal panel weighting and target/family-aware uncertainty. Do not let a large mutation family or duplicated model rows dominate the evaluation. Report excluded/failed panels alongside successful ones.

Mutation panels test sensitivity to affinity changes among related variants. They do not by themselves establish performance for unrelated de novo binder scaffolds or calibrated absolute affinity prediction.

### 7.4 Keep computation bounded

A future first pilot could use **3–5 eligible panels with 6–10 candidates each**, selected before model outcomes are known. This is a planning proposal, subject to eligibility and measured runtime.

With five predictions per ordered input including self-pairs, that design spans roughly **540–2,500 predictions per engine**. A single 100-candidate panel would require **50,000 predictions per engine**. These counts demonstrate why the larger benchmark should not become an implicit obligation of this $0 first release.

Before starting that follow-up, produce an eligibility manifest, panel counts, explicit prediction counts, and a compute estimate using measured timings. The current plan does not authorize paid compute.

## 8. Validation gates, release sequence, and completion

### 8.1 Required gates

| Gate | Required evidence | Effect on release |
|---|---|---|
| Historical scoring | Archived CSV verified and reproduced through the normal entry point; justified correctness changes documented. | Required for the scientific regression baseline. |
| Package and CPU path | Candidate wheel works outside the checkout, bundled examples load, targeted checks and headless quantification pass. | Required. |
| Shared runner | Python/command/notebook use shared functions; resume and failure semantics pass focused tests. | Required. |
| Public guided modeling path | Candidate notebook reaches rankings and heatmap in genuine free Colab using the proposed default or supported fallback. | Required for the corresponding public default. |
| ESMFold2 recommendation | Frozen full six-nanobody gate and end-to-end/runtime requirements pass. | Required to promote ESMFold2; otherwise fallback/experimental/blocked policy applies. |
| Publication | Released wheel installed afresh; package/version/artifact provenance and public notebook links verified. | Required to call the public release complete. |
| Larger affinity benchmark | Recommendation and follow-up protocol documented. | No large-dataset folding gate for this release. |

### 8.2 Prepare the release candidate

Check current PyPI history and choose an unused version. The initial proposal remains `1.3.2`; do not reuse an already published or yanked version. [PyPI history](https://pypi.org/project/appraise/#history)

Build wheel and source distribution from a clean candidate checkout, run `twine check`, inspect packaged code/data, and install the wheel outside the checkout. Verify the scientific baseline, saved-structure example, runner, and notebook/package alignment. Use the exact candidate artifact or revision in Colab.

Create or update CPU CI in the canonical repository. Link real Colab logs and the scientific comparison from the PR; distinguish unrun checks from passing ones. Avoid committing weights or large generated structures.

### 8.3 Merge and publish in a working order

1. Before the first merge, decide whether notebook activation must be separate. If it imports candidate-only package symbols, prepare a tested activation change while keeping the current public notebook on a working released package. Establish the actual publisher trigger and artifact handoff.
2. Complete required checks and merge the package/runner changes through the canonical repository's normal protections, preserving the working public notebook during the publication interval.
3. Build release artifacts from the verified merged revision, verify code/data contents against the tested candidate, install/test the actual release artifacts, and record hashes. The publisher must upload these artifacts or validate its rebuilt artifacts before upload.
4. Invoke the configured publication route in its required order. For a tag-triggered publisher, create the matching tag on the verified revision to start the workflow; for a manual workflow, follow its declared trigger. Verify the PyPI publication before announcing the completed GitHub release. Reconcile already-published files/tags before retrying an interruption.
5. Merge the already-tested public notebook activation pointing to the exact published package. Verify that README and notebook badges launch the canonical `xz-ding/APPRAISE` notebook.
6. Keep synchronization of the laboratory fork outside the critical publication path.

Prepare a concrete Trusted Publisher/workflow configuration before any owner-only PyPI setup handoff. Follow the [PyPI publishing documentation](https://docs.pypi.org/trusted-publishers/).

### 8.4 Verify as a new user

Read back the published version and hashes. Install from PyPI in a fresh environment and verify the saved-structure and runner paths. Open the default branch's actual public Colab link and run the guided example from installation to results on a qualified free configuration.

This is the public usability check. Candidate-wheel tests do not establish that public links or publication resolve correctly.

### 8.5 Final completion record

Report:

- canonical merged PR(s), commit, release tag, and PyPI version;
- public notebook URL and verified package/artifact hashes;
- CPU and actual-Colab validation evidence;
- exact default-engine policy and measured operating range;
- six-nanobody ranks, correlation, coverage, and qualification outcome;
- support labels for retained advanced engines;
- durable checkpoint and exact next action for any external blocker;
- the separately scoped SKEMPI follow-up recommendation.

If ESMFold2 misses the gate, report the measured result and ship the qualified fallback when its release gates pass. If ESMFold2 cannot be exercised, report unresolved qualification explicitly. If publication cannot proceed, report the concrete publication blocker without claiming that the public upgrade is finished.

