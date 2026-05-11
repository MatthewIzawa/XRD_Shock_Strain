---
title: Reference-phase sensitivity test for W-A (JAC revision Table S1) — design
date: 2026-05-11
status: design
scope: Llunr/refphase_sensitivity/ (analysis) → JAC Table S1 + §4.1 prose (paste-into-docx)
related: addresses Kleppe reviewer comment 4
author: Matthew R. M. Izawa
---

# Reference-phase sensitivity test for W-A — design

## 1. Background

Reviewer comment (Kleppe):

> "You may want to briefly comment on the sensitivity of the W-A
> results to reference phase selection (e.g. anorthite choice),
> especially for compositionally variable samples."

The JAC manuscript reports plagioclase column lengths derived from
Warren-Averbach analysis using pure anorthite (CaAl₂Si₂O₈) as the
reference phase. The Plag/Pyx column-length ratio crossover at the
maskelynite-forming shock threshold (Figure 2) is the headline
methodological result. A reviewer can reasonably ask: does this
result hold when the reference-phase composition mismatches the
actual plagioclase in the sample?

This spec defines a focused sensitivity test comparing two endmember
plagioclase references — pure albite (NaAlSi₃O₈) and pure anorthite —
for three eucrites spanning the shock series, sufficient to retire
the reviewer comment without expanding manuscript scope.

## 2. Goal

Produce one supplementary table (Table S1) and ~3 sentences of prose
for §4.1, quantifying the W-A method's sensitivity to plagioclase
reference-phase choice.

## 3. Scope

**In scope:**

- Three samples: Tirhert (S1), NWA 1836 (S3), JaH 626 (S6) —
  eucrite shock series, all measured at I11.
- Two plagioclase references: anorthite (existing CIF) and low
  albite (Harlow & Brown 1980, AMCSD 0000797, already placed in
  `Llunr/CIFs/Albite__0000797.cif`).
- One supplementary table comparing D_med (Å) and family count for
  each (sample, reference) pair.
- One paragraph (~3 sentences) of §4.1 prose, contingent on outcome.
- One reproducible analysis script in a new standalone
  `Llunr/refphase_sensitivity/` directory (parallel to
  `Llunr/comparison_v040/`), keeping `Paper1_JAC/` untouched.

**Out of scope (deferred or rejected):**

- Pyroxene reference sensitivity — separate question; diopside and
  enstatite CIFs are on hand for future work but not exercised here.
- Intermediate plagioclase compositions (bytownite, oligoclase) —
  endmember comparison is the sharpest reviewer-response.
- Sensitivity across the full 29-sample survey — three eucrites
  span the shock-relevant compositional range, and the headline
  result lives in the eucrite series.
- Changes to the `xrd_profile` package itself — uses read-only API.
- Changes to existing JAC pipeline scripts (`compile_survey.py`,
  `figures_jac.py`, etc.) or the `survey_results_29samples.csv`.

## 4. Inputs (all pre-existing)

- I11 raw XRD data, .xye format, for the three samples (in
  `Desktop/111 Backup.../IPM/2018/ee17803-1/processing/`,
  λ = 0.826517 Å). Read-only.
- Anorthite reference:
  `Llunr/xrd_profile/examples/cifs/Anorthite.cif`. Read-only.
- Albite reference:
  `Llunr/CIFs/Albite__0000797.cif` (Harlow & Brown 1980, low
  albite, Amelia Virginia type locality, near-endmember
  Na₀.₉₈₆Al₁.₀₀₅Si₂.₉₉₅O₈, C-1 triclinic). Read-only.
- `xrd_profile` v0.4.1, installed in the Anaconda base env.

## 5. Approach

### 5.1 Alternatives considered (and rejected)

**Alt A — extend `compile_survey.py`.** Add an albite-reference
column to the master survey CSV. Rejected: mixes a sensitivity-test
artifact into the frozen 29-sample survey, violates the freeze
decision recorded in `feedback_no_overwrite_submission_figures.md`
and `project_shock_pipeline_decisions.md`.

**Alt B — add a `sensitivity_compare()` utility to xrd_profile.**
PRO: would be a reusable package feature. Rejected: this is a
paper-specific one-off; the package already exposes `Phase.from_cif`
and `XRDProfile.run_all`, which are sufficient. Adding sensitivity
machinery to v0.4.1 also pollutes the release boundary.

**Alt C (chosen) — standalone `Llunr/refphase_sensitivity/`
directory.** Single self-contained
`tableS1_reference_sensitivity.py` consuming the xrd_profile
public API read-only, writing only into its own directory.
Parallel to the existing `Llunr/comparison_v040/` pattern. No
package change. No `Paper1_JAC/` contamination. Final deliverables
(Table S1 markdown and §4.1 prose) reach the manuscript through
paste-into-docx by the author, never touching `Paper1_JAC/` files
directly.

### 5.2 Per-sample logic

For each (sample, reference) of the six combinations:

1. Load the I11 .xye file via `XRDProfile.from_file(...,
   wavelength=0.826517)`.
2. Load the reference CIF via `Phase.from_cif(...)`.
3. Run W-A via `XRDProfile.guided_warren_averbach(phase=...)` or
   equivalent (preferred entry point of v0.4.1 — implementer to
   confirm).
4. Capture from the result object:
   - Median column length `D_med` (Å)
   - Number of valid harmonic families passing the quality filter
   - Number of detected peaks (diagnostic, retained in CSV but not
     in the published table)
   - Per-family D values (diagnostic, retained in CSV only)

## 6. Outputs

### 6.1 CSV (full record)

`Llunr/refphase_sensitivity/tableS1_reference_sensitivity.csv` —
six rows (3 samples × 2 references), columns: `sample`, `shock`,
`reference_phase`, `cif_path`, `n_peaks_detected`, `n_families`,
`D_med_A`, `notes`. A per-peak diagnostic file is written
alongside (see §7).

### 6.2 Supplementary table (published form)

| Sample (shock) | D_med, anorthite (Å) | D_med, albite (Å) | Δ (%) | Families, anorthite | Families, albite |
|---|---:|---:|---:|---:|---:|
| Tirhert (S1) | _to fill_ | _to fill_ | _to fill_ | _to fill_ | _to fill_ |
| NWA 1836 (S3) | _to fill_ | _to fill_ | _to fill_ | _to fill_ | _to fill_ |
| JaH 626 (S6) | _to fill_ | _to fill_ | _to fill_ | _to fill_ | _to fill_ |

Δ% = 100 × |D_med(albite) − D_med(anorthite)| / D_med(anorthite).

### 6.3 Prose for §4.1 (contingent on outcome)

Two templates, one chosen after the script runs.

**Low-sensitivity case (mean Δ < ~10 %):**

> "We tested the sensitivity of the W-A column lengths to
> plagioclase reference-phase choice by repeating the analysis for
> three eucrites (Tirhert S1, NWA 1836 S3, JaH 626 S6) using pure
> albite (NaAlSi₃O₈; Harlow & Brown, 1980) in place of pure
> anorthite (CaAl₂Si₂O₈). Median column lengths agreed within X %
> and family counts within ±Y across the compositional endpoints
> of the plagioclase solid solution (Table S1), confirming that
> the Plag/Pyx ratio reported above is robust to the precise
> reference CIF used."

**Moderate-/high-sensitivity case (mean Δ ≥ ~10 %):**

> "We tested the sensitivity of the W-A column lengths to
> plagioclase reference-phase choice by repeating the analysis for
> three eucrites (Tirhert S1, NWA 1836 S3, JaH 626 S6) using pure
> albite (NaAlSi₃O₈; Harlow & Brown, 1980) in place of pure
> anorthite (CaAl₂Si₂O₈). Median column lengths varied by X % and
> family counts by ±Y across the plagioclase compositional
> endpoints (Table S1). For compositionally variable samples,
> composition-matched reference CIFs are therefore recommended when
> absolute column lengths are quoted; the relative Plag/Pyx ratio
> derived from a single self-consistent reference remains an
> internally valid shock indicator."

## 7. Risks and mitigations

- **Peak-detection truncation.** If the reference-guided peak
  detection's d-spacing tolerance window is narrower than the
  albite-vs-anorthite peak-position offsets for this sample
  composition, the albite run will detect *fewer peaks*, and the
  family count drop will partly reflect "missed peaks" rather than
  "different reference." Mitigation (per the user, accepted): the
  implementation will export a per-peak diagnostic table for each
  (sample, reference) pair — peak position, intensity, family
  assignment, and the d-spacing offset between the detected peak
  and the nearest reference-predicted peak. The implementer
  reviews these diagnostics to confirm detected peaks correspond
  to actual sample peak maxima under each reference, and flags
  any wrongly-attributed peaks in the §4.1 prose. The CSV also
  captures `n_peaks_detected` separately from `n_families`; if
  they diverge qualitatively, the prose distinguishes the two
  effects.

- **Wrong API entry point.** v0.4.1 exposes several W-A entry
  points (`run_all`, `guided_warren_averbach`, etc.) Implementer
  to confirm the canonical entry point matches what
  `compile_survey.py` uses, so the comparison is apples-to-apples.

- **Symmetry / lattice-parameter compatibility.** Both phases are
  triclinic feldspars (low albite C-1, anorthite P-1 with C-1
  pseudo-cell). Compatible. Lattice-parameter offsets are modest
  (a, b, c within ~5 %; angles within a few degrees). Pymatgen will
  generate distinct reference peak lists, which is the whole point.

## 8. Acceptance criteria

- Script runs end-to-end on all six combinations without errors.
- CSV populated with six rows and the eight columns above.
- Markdown table populated from the CSV, with Δ% computed.
- One of the two prose templates adapted with concrete numbers
  and presented for paste-into-docx.
- All edits respect the JAC freeze: no changes to the 29-sample
  survey, no changes to existing figures, no changes to
  `xrd_profile` source, no new files added to `Paper1_JAC/`.
- Per-peak diagnostic file written alongside the main CSV;
  implementer reviews it before drafting prose.
- New script saved to `Llunr/refphase_sensitivity/` (outside any
  git repo, so "committed" here means saved to disk and
  reproducible from the spec).

## 9. Deliverables

1. `Llunr/refphase_sensitivity/` — new standalone directory,
   parallel to `Llunr/comparison_v040/`. Created if absent.
2. `Llunr/refphase_sensitivity/tableS1_reference_sensitivity.py` —
   single self-contained analysis script (~150 lines).
3. `Llunr/refphase_sensitivity/tableS1_reference_sensitivity.csv`
   — script output (six rows; main results table).
4. `Llunr/refphase_sensitivity/per_peak_diagnostics.csv` —
   diagnostic file (per-peak position, intensity, family
   assignment, d-spacing offset from reference prediction) used
   for the §7 mitigation review.
5. Markdown table text (Table S1) — printed to stdout and
   presented in chat for paste into the .docx supplementary.
6. §4.1 prose text — one of two templates filled in with the
   numerical result, presented in chat for paste into the .docx
   body.

## 10. Adjacent revision items (separately addressable)

- **Kleppe comment 2** (sharpen Introduction to emphasize the
  methods-paper / validation-demo framing): small prose edit, no
  scripts. Can be drafted in the same revision pass via a separate
  text deliverable in chat.
- **Kleppe comment 3** (W-H multi-phase limitations as a strength):
  praise, no action required.
- **Kleppe comment 1** (compelling Plag/Pyx take-home): praise,
  no action required.
