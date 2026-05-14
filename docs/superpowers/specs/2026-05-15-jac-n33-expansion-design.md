---
title: JAC N=33 expansion — capillary subtraction + 4 new samples + 2 lab→synchrotron upgrades — design
date: 2026-05-15
status: design
scope: Llunr/n33_expansion/ (new analysis) + Paper1_JAC/ (manuscript, figures, tables, cover letter)
related: post-co-author-circulation DLS data arrival (24 Apr 2026); strengthens §3.3 Plag/Pyx reversal; retires §4.3 capillary "would enable" sentence
author: Matthew R. M. Izawa
---

# JAC N=33 expansion — design

## 1. Background

On 13 May 2026 a Diamond Light Source I11 MAC dataset from the 24 Apr 2026
beamtime (proposal `Izawa_achondrites_20260424`) arrived at
`Llunr/HED_XRD_Shock/New_DLS_20260513/`. It contains:

- 6 sample scans, ready-to-process 2θ-counts-error text format
  (1437211 NWA 7465, 1437212 NWA 5478, 1437213 NWA 5751, 1437214 NWA 7831,
  1437215 NWA 6693, 1437218 Talampaya).
- 1 blank-capillary background scan (1437220).
- 3 entries in the scan log have no corresponding `.dat` file (1437216
  Brachina, 1437217 NWA 6013, 1437219 Tatahouine — failed or not yet
  processed; no impact on this spec).

The current JAC manuscript (`Paper1_JAC/manuscript_paper1_JAC.md`) is in
post-co-author-circulation state with the PJAM-MI version finalised on
9 May 2026, the abstract rewritten to the asymmetric-amorphisation
framing, and the cover letter rewritten with Greek α and no
corporate-speak. Submission was imminent.

Two of the new data items directly affect the manuscript's standing claims:

1. **The blank capillary measurement (1437220).** §4.3 of the manuscript
   currently states:

   > "...the borosilicate capillary used for synchrotron measurements
   > contributes an amorphous signal that overlaps with shock-produced
   > glass. A blank capillary measurement at the same beam conditions
   > would enable subtraction of this contribution."

   This is now a present-tense capability. Shipping the manuscript with
   the "would enable" sentence intact reads as dated by the time
   reviewers see it.

2. **NWA 5478 (shocked polymict eucrite) and Talampaya upgrade.** The
   §3.3 Plag/Pyx reversal currently rests on two high-shock eucrites
   (Millbillillie S5, JaH 626 S6). NWA 5478 provides a third
   independent high-shock data point. Talampaya, currently in the
   survey as a Lab Co cumulate eucrite (S2), now has a synchrotron
   measurement that adds another low-shock eucrite to the reversal's
   low-shock anchor (currently Tirhert S1 and Bereba S2-S3).

The user has approved a pre-submission delay to integrate the new data.
This spec defines the integration.

## 2. Goal

Produce a v2 JAC manuscript with `N = 33` samples spanning 8 parent-body
groups, capillary-subtracted PDFs across all 18+ synchrotron samples,
and a strengthened §3.3 Plag/Pyx reversal claim anchored by 3 high-shock
eucrites instead of 2. Submit v2 to co-authors for a single approval
round before JAC submission.

## 3. Scope

### 3.1 In scope

- Process the 7 new `.dat` files and re-process the 18 existing I11
  patterns with capillary subtraction.
- Run the `xrd_profile` v0.3.x analysis pipeline on the 6 new sample
  scans (4 new samples + 2 lab→synchrotron upgrades) to extract W-A
  column lengths and PDF first-peak FWHMs.
- Re-run PDF first-peak fits on the 18 existing I11 samples with
  capillary-subtracted intensities.
- Regenerate the `survey_results_*.csv` as `survey_results_33samples.csv`.
- Regenerate every data-bearing figure (Figs 1, 2, 3, 4, 5, 6, 7) into
  `Paper1_JAC/paper1_figures/`.
- Update the manuscript text in `Paper1_JAC/manuscript_paper1_JAC.md`
  for the new sample counts, group count (8), capillary subtraction,
  and strengthened Plag/Pyx claim.
- Update Table 1 (29 rows → 33 rows; 18+ PDF FWHM values; 2 row replacements),
  Table 2 (3 new eucrite rows including NWA 5478), and add a new
  Table S2 (cross-instrument comparison for NWA 5751 and Talampaya).
- Update cover letter, suggested reviewers, acknowledgments, data
  availability statement.
- Generate v2 .docx with track changes against the PJAM-MI v1.
- Single round of co-author re-circulation.

### 3.2 Out of scope

- **No `xrd_profile` package code changes.** Capillary subtraction is
  implemented as a standalone analysis script in `Llunr/n33_expansion/`
  that consumes the existing `xrd_profile` API. The package stays at
  v0.3.x for the manuscript-of-record; the v0.4.0/v0.5.0/v1.0.0 work
  planned elsewhere is independent.
- **No changes to `Paper1_JAC/submission/figures/`** (Inkscape
  hand-edited; per `~/.claude/projects/.../MEMORY.md`). All Python
  figure regeneration outputs into `Paper1_JAC/paper1_figures/`; the
  user manually decides which to re-import via Inkscape.
- **No structural figure additions to the main body** beyond what
  already exists. Optional new figures (capillary before/after
  demonstration) live in supplementary.
- **No new co-author additions** even though the new beamtime may
  warrant one. Authorship locked to the PJAM-MI 11-author list.
- **No re-analysis of W-A column lengths** on the existing 18
  synchrotron samples. Capillary contributes to amorphous diffuse,
  not crystalline Bragg peak profiles; W-A is unaffected.
- **NWA 6013 and Tatahouine remain Lab Co rows.** Their Apr 2026 I11
  attempts failed; they are not part of the upgrade pair.
- **Re-analysis with newer `xrd_profile` releases.** This spec uses
  the v0.3.x API exclusively to keep the manuscript reproducible at
  one package tag.

### 3.3 Decided during brainstorm

| Decision | Choice |
|---|---|
| Scope of new-data inclusion | Full inclusion (N=33) |
| Timing | Pre-submission delay, 5-phase validate-first |
| Architecture | Option B: validate methods first, then expand |
| Parent-body groups | 8 (NWA 6693 is a new ungrouped achondrite group, paired with NWA 6704; oxidised assemblage with Fe-rich olivine, Fe-rich pyroxene, awaruite trace; low shock) |
| Capillary subtraction scaling | High-Q tail fit (Q ≈ 12-14 Å⁻¹, sample-dependent) |
| Capillary mention in abstract | No — this is a software paper; capillary subtraction lives in §2.5 / §4.3 only |
| Fig 5 cross-instrument bonus | Yes — plot NWA 5751 and Talampaya on both Lab Co and I11 panels with connecting lines |
| Fig 7 before/after | Likely supplementary, not main text |
| Cross-instrument supplementary table | Yes — new Table S2, Lab Co vs I11 W-A column lengths for NWA 5751 and Talampaya |

## 4. Inputs

### 4.1 New (this batch)

- `Llunr/HED_XRD_Shock/New_DLS_20260513/1437211-mac-summed.dat` — NWA 7465
- `…/1437212-mac-summed.dat` — NWA 5478
- `…/1437213-mac-summed.dat` — NWA 5751 (upgrade)
- `…/1437214-mac-summed.dat` — NWA 7831
- `…/1437215-mac-summed.dat` — NWA 6693
- `…/1437218-mac-summed.dat` — Talampaya (upgrade)
- `…/1437220-mac-summed.dat` — blank capillary

All 7 files: I11 MAC, 1800 s scans, 24 Apr 2026, 150 020 2θ points,
`Wavelength = Not Set` in header (must be confirmed; presumed 0.8265 Å
matching previous I11 beamtime).

### 4.2 Pre-existing

- `Paper1_JAC/manuscript_paper1_JAC.md` (PJAM-MI v1, 11 authors,
  N=29, 7 groups).
- `Paper1_JAC/survey_results_29samples.csv` (current canonical
  per-sample analysis output).
- `Paper1_JAC/table1_sample_inventory.md` (current Table 1 source).
- `Paper1_JAC/submission/cover_letter.md`,
  `…/suggested_reviewers.md`, `…/build_auxiliary_docx.py`.
- `Paper1_JAC/paper1_figures/` (current PNG outputs from regeneration
  scripts).
- `Paper1_JAC/submission/figures/` — **off-limits, Inkscape hand-edited**.
- `xrd_profile` v0.3.x, installed in the Anaconda base env.

## 5. Phase 1 — Methods validation (1-2 days)

**Goal.** Prove the capillary subtraction does not undermine the PDF
FWHM monotonic-with-shock trend before any further commitment.

**Steps.**

1. **Confirm beamtime wavelength.** Check the 24 Apr 2026 I11 MAC
   proposal documentation / DLS log book. If λ ≠ 0.8265 Å the
   downstream Q-space numerics shift; stop and re-scope.
2. **Reduce the capillary scan to `.xy` form** matching existing
   I11 file convention. Sanity check: broad amorphous hump only,
   no crystalline Bragg peaks.
3. **Implement capillary subtraction with high-Q-tail scaling**:
   - Fit a scale factor `s` minimising
     `Σ_(Q ∈ [12, 14]) (I_sample(Q) − s · I_cap(Q))²`
     over a sample-dependent high-Q window where the sample's
     coherent Bragg contribution is negligible.
   - Subtract: `I_corrected(Q) = I_sample(Q) − s · I_cap(Q)`.
4. **Pick 3 test eucrites spanning shock**: Tirhert (S1), NWA 6477 (S3),
   JaH 626 (S6).
5. **Run `xrd_profile` PDF on capillary-subtracted intensities** for
   each test sample. Extract first-peak FWHM via the existing
   Gaussian-fit module.
6. **Produce one diagnostic figure** (4 panels):
   - (a) Capillary spectrum I(Q) alone.
   - (b) Sample I(Q) before and after capillary subtraction for the
     3 test samples, low-Q region emphasised.
   - (c) PDF G(r) before vs after for the 3 test samples,
     r = 0-10 Å.
   - (d) Full FWHM-vs-shock trend across all 18 existing I11 samples
     with the 3 test samples updated (open vs filled markers).
7. **Decision gate.** Author reviews diagnostic figure.
   - GO: monotonic broadening trend preserved or strengthened;
     FWHM rank-order with shock stage intact.
   - PAUSE: trend disrupted. Decide whether to (i) refine the
     scaling method, (ii) reinterpret the FWHM trend, or (iii) ship
     capillary subtraction as a supplementary methods improvement
     only.

**Deliverable.** `Llunr/n33_expansion/phase1_validation/diagnostic.png`
+ a 3-paragraph note recording the GO/PAUSE decision and the chosen
scaling-method parameters.

**Abort condition.** If the trend is disrupted in a way that cannot
be resolved by refining the scaling method, fall back to Option A
("Capillary only methods improvement") or escalate to the user for
a re-scoping discussion.

## 6. Phase 2 — Data processing and survey expansion (~1 week)

**Goal.** Produce `survey_results_33samples.csv` with capillary-subtracted
PDFs and the 6 new analyses.

**Steps.**

1. **Verify pre-expansion inventory against the canonical CSV**
   (`Paper1_JAC/survey_results_29samples.csv`): confirm 7 Lab Cu +
   18 I11 + 4 Lab Co = 29, and the I11 breakdown of 12 HED + 6 non-HED.
   Any discrepancy with this spec's numbers gets resolved before
   building the new CSV. Then reduce all 7 new `.dat` files to `.xy`
   matching the existing I11 file convention. Place under
   `Llunr/n33_expansion/data/i11_2026Apr_reduced/`.
2. **Apply the Phase-1 validated capillary subtraction** to:
   - the 6 new sample scans
   - the 18 existing I11 sample patterns
   Cache the capillary-subtracted intensities as `.xy` files in
   `Llunr/n33_expansion/data/i11_corrected/` (existing samples
   inherit the same filename pattern with `_cap` suffix or similar).
3. **Run `xrd_profile.XRDProfile.run_all` on the 6 new analyses**.
   For each: W-A column-length statistics (per-family ⟨L⟩_area,
   sample-level median Bertaut intercept, harmonic family count, etc.),
   PDF first-peak Gaussian fit (FWHM, position), Williamson-Hall
   (recorded for completeness). Output per-sample JSON or CSV results
   into `Llunr/n33_expansion/results/per_sample/`.
4. **Re-run the PDF first-peak fit on the existing 18 I11 samples**
   with capillary-subtracted intensities. PDF only; do not re-run W-A
   (capillary contributes to amorphous diffuse, not crystalline Bragg
   peak profiles).
5. **Sample metadata research** for the 4 new samples (parallel task,
   can run concurrently with steps 1-4):
   - NWA 7465 — monomict eucrite: lit shock stage, modal mineralogy,
     primary reference (Meteoritical Bulletin entry).
   - NWA 5478 — shocked polymict eucrite: lit shock stage (likely S5
     or S6), modal mineralogy, reference.
   - NWA 7831 — diogenite: lit shock stage, opx vs olv ratio, reference.
   - NWA 6693 — CC-associated impact melt / NWA 6704-paired,
     ungrouped achondrite: low shock, oxidised assemblage (Fe-rich
     olivine, Fe-rich pyroxene, awaruite Ni₃Fe trace), primary
     reference (Warren et al. is a candidate, e.g., the 6704 paper).
   - `[AUTHOR INPUT NEEDED]` placeholders acceptable for v2 co-author
     circulation; verify references before final JAC submission.
6. **Aggregate into `survey_results_33samples.csv`**:
   - 4 new rows (4 new samples, all I11)
   - 2 rows replaced (NWA 5751, Talampaya — Lab Co rows out, I11 rows in)
   - 18 rows updated for PDF FWHM column (existing I11 samples)
   - 13 rows untouched (7 lunar Lab Cu + 2 remaining Lab Co HEDs +
     4 non-HED I11 samples already counted in the 18)
   - **Inventory check:** 7 Lab Cu + 17 HED I11 (12 existing + 2
     upgrades + 3 new: NWA 7465, 5478, 7831) + 7 non-HED I11
     (6 existing + NWA 6693 new) + 2 HED Lab Co (NWA 6013,
     Tatahouine, I11 attempts failed) = **33 ✓**.
7. **Preserve cross-instrument data**: keep the existing Lab Co rows
   for NWA 5751 and Talampaya in a parallel CSV
   `Llunr/n33_expansion/results/cross_instrument_pairs.csv` for
   downstream Table S2 generation. These rows are *not* in the main
   Table 1 but are real measurements.
8. **Phase 2 sanity-check gate.** Inspect the 6 new sample analyses:
   - Does NWA 5478 (shocked polymict eucrite) have Plag/Pyx column-length
     ratio < 1, consistent with the reversal prediction?
   - Does Talampaya (cumulate eucrite S2) have synchrotron W-A column
     length larger than its Lab Co value (310 Å), and Plag/Pyx > 1?
   - Does NWA 7831 (diogenite) sit reasonably alongside Tatahouine,
     NWA 6013, NWA 2968?
   - Does NWA 6693 sit at the low-shock end as expected for an
     unshocked CC-associated impact melt?
   - Author reviews. If anomalous, understand why before regenerating
     figures.

**Deliverable.** `survey_results_33samples.csv` + a one-page numbers-check
table comparing the 6 new analyses against predictions.

**Abort condition.** If the sanity check reveals an anomaly that
implies a measurement or processing error, halt and resolve before
Phase 3.

## 7. Phase 3 — Figure regeneration (~2-3 days)

**Goal.** Every data-bearing figure reflects N=33 with capillary-subtracted
PDF values. All regeneration outputs land in `Paper1_JAC/paper1_figures/`;
`Paper1_JAC/submission/figures/` stays untouched.

**Per-figure changes.**

- **Fig 1.** If panels are "one per parent-body group," add an 8th
  panel for NWA 6693. If panels are scatter-plot, re-render with new
  symbol mapping and updated counts. Caption sample count 29→33,
  group count 7→8.
- **Fig 2 (anpig crossover, Plag/Pyx vs shock).** Add 3 eucrites:
  NWA 7465 (mid), NWA 5478 (high — new third high-shock anchor),
  Talampaya (S2 low — synchrotron upgrade). Plag/Pyx reversal gets
  a real third high-shock data point.
- **Fig 3 (PDF FWHM vs W-A column length, synchrotron only).** All 18
  existing synchrotron points get new FWHM (capillary-subtracted);
  +6 new points. Total 24 synchrotron points. New symbol/legend
  entry for NWA 6693 (8th parent-body group).
- **Fig 4 (W-H reliability classification).** Add the 6 new samples to
  the classification for completeness. Caption count update.
- **Fig 5 (families vs Q_max, cross-instrument).**
  - Panel (a) Lab Cu: 7 samples → 7 (unchanged).
  - Panel (b) Lab Co: 4 samples → 2 (NWA 6013, Tatahouine).
  - Panel (c) I11: 18 samples → 24.
  - **Cross-instrument bonus (approved):** plot NWA 5751 and
    Talampaya on both panels (b) and (c), with connecting lines
    between the same-specimen pair, demonstrating Q_max effect on
    the same physical sample.
- **Fig 6 (W-A median column length vs first PDF peak position).** All
  I11 PDF positions and FWHMs update with capillary subtraction;
  +6 new points. Confirm Tagish Lake serpentine outlier and NWA 6711
  anomalous first-shell still separate cleanly.
- **Fig 7 (PDF examples / G(r) curves).** All curves shown become
  capillary-subtracted. Optionally add NWA 5478 and/or NWA 6693 panels.
- **Fig S1 (column-length schematic).** Unchanged.
- **Fig S2 (new, supplementary).** Capillary before/after
  demonstration: one sample's G(r) before vs after subtraction +
  capillary spectrum itself. Demonstrates why the new measurement
  matters.

**Phase 3 sanity-check gate.**

- Fig 2: Plag/Pyx reversal still visible, now anchored by 3 high-shock
  eucrites?
- Fig 3 / Fig 6: PDF FWHM monotonic-with-shock trend still visible,
  ideally cleaner?
- Fig 5: cross-instrument pairs for NWA 5751 / Talampaya look as
  expected (more families at higher Q_max)?

**Deliverable.** Updated PNGs in `Paper1_JAC/paper1_figures/` + a
side-by-side comparison sheet (old vs new) + a punch list of which
figures need Inkscape rework before they re-enter `submission/figures/`.

## 8. Phase 4 — Text and metadata updates (~3-4 days)

**Goal.** Every text reference to sample counts, instruments, group
counts, and PDF backgrounds is coherent for N=33, capillary subtraction,
and 8 parent-body groups.

### 8.1 Sample-metadata research (parallel task)

For each of the 4 new samples: Meteoritical Bulletin entry, lit shock
stage, modal mineralogy, primary reference. `[AUTHOR INPUT NEEDED]`
placeholders acceptable for v2 co-author circulation; verify before
final JAC submission.

### 8.2 Manuscript prose updates, section by section

| Section | Change |
|---|---|
| Abstract | "29-meteorite suite spanning seven parent-body groups" → "33-meteorite suite spanning eight parent-body groups". Capillary subtraction *not* added (per decision). |
| §1 Intro | Search-and-replace any hardcoded "29" or "seven". |
| §2.5 PDF methods | Add a short paragraph: blank-capillary subtraction with high-Q-tail scaling, applied before Chebyshev background, Lorch modification, and Gaussian first-peak fit. |
| §2.6 Samples and data collection | Sample count 29→33. Add 4 new samples to the list with instrument (I11 MAC), date (Apr 2026), exposure. Mention upgrades (Talampaya, NWA 5751 measured on both instruments). Mention the blank capillary measurement. |
| §3.1–§3.x Results | Number updates throughout. §3.3 (Plag/Pyx reversal) gets the most attention: NWA 5478 as a third high-shock anchor; Talampaya joins the low-shock side (with cross-instrument note); NWA 7465 fills the middle. |
| §4.1 (albite/anorthite sensitivity) | Unchanged. |
| §4.2 (asymmetric-amorphisation rescue) | Short addition: NWA 5478 strengthens the third high-shock anchor. |
| §4.3 (PDF limitations) | **The "would enable" sentence about the capillary is replaced** by a sentence describing the subtraction now in place, with results referenced in §3. Keep the lab-Q_max-vs-synchrotron caveat. |
| §4.4 Software design | Unchanged. |
| §5 Conclusions | Number updates. Asymmetric-amorphisation conclusion gets the third high-shock anchor. Modify or add a conclusion on capillary subtraction. |

### 8.3 Tables and references

- **Table 1.** 4 new rows; 2 rows replaced (NWA 5751, Talampaya → I11
  row); 18 rows updated (PDF FWHM column). Caption sample count.
- **Table 2** (Plag/Pyx ratios for eucrites). 3 new rows
  (NWA 7465, NWA 5478, Talampaya). NWA 5478 is the key new row.
- **Table S1** (albite/anorthite sensitivity). Unchanged.
- **Table S2** (new). Cross-instrument comparison: Lab Co vs
  synchrotron W-A column lengths for NWA 5751 and Talampaya.
- **References.** Add Meteoritical Bulletin / classification refs for
  the 4 new samples. Add a reference for capillary subtraction in
  synchrotron PDF practice if appropriate (Billinge group or
  Wright et al.).

### 8.4 Submission documents

- **Cover letter.** "29-meteorite suite spanning seven parent-body
  groups" → updated for N=33 and 8 groups. Capillary subtraction
  not mentioned (consistent with abstract decision).
- **Suggested reviewers.** Profile B (PDF methodology) text remains
  accurate. No structural changes; re-circulate for completeness.
- **Acknowledgments.** Add the Apr 2026 Diamond Light Source
  proposal/beamtime number.
- **Data availability statement.** Add the new `.dat` files to the
  archived dataset; bump Zenodo DOI version if applicable.

**Deliverable.** Updated `manuscript_paper1_JAC.md` + updated cover
letter + updated suggested reviewers + updated data availability +
sample-metadata punch list.

## 9. Phase 5 — Co-author re-circulation (~2-3 weeks calendar)

**Goal.** Present v2 to the 11 co-authors as a single clean update with
track changes against v1 (PJAM-MI) and a cover note that orients them
in 60 seconds.

**Steps.**

1. **Generate v2 .docx with track changes** against v1. Run
   `Paper1_JAC/build_docx.py` on the updated Markdown to produce
   `manuscript_paper1_JAC_v2.docx`. Use Word's *Compare* (or
   pandoc track-changes diff) against
   `manuscript_paper1_JAC_submission_EAC may 9 2026-PJAM-MI.docx`.
   Apply the standard .docx hardening per `~/.claude/CLAUDE.md`
   §8 (no em-dashes in prose, no coloured headings, black-only,
   Unicode super/subscripts → `vertAlign`, captions/footnotes 10pt).
   Pin v1 as `..._v1_approved.docx`.

2. **Prepare co-author cover note** (one page).
   - What's new in v2:
     - 4 new I11 samples (NWA 7465, NWA 5478, NWA 7831, NWA 6693)
     - 2 lab→synchrotron upgrades (NWA 5751, Talampaya)
     - Blank capillary subtraction (validation in Phase 1, integrated
       in §2.5 and §4.3)
     - All data-bearing figures regenerated
     - §4.3 "would enable" sentence retired
   - What strengthened: Plag/Pyx reversal anchored by 3 high-shock
     eucrites; PDF FWHM trend now capillary-corrected; cross-instrument
     caveat now backed by same-specimen measurements (NWA 5751,
     Talampaya).
   - What's unchanged: authorship and order; suggested reviewers;
     asymmetric-amorphisation interpretation in §4.2; `xrd_profile`
     methodology.
   - Outstanding `[AUTHOR INPUT NEEDED]` items, by section, with
     deadline.

3. **Bundle for distribution.**
   - `manuscript_paper1_JAC_v2_trackchanges.docx`
   - `manuscript_paper1_JAC_v2_clean.docx`
   - `cover_letter_v2.docx`
   - `suggested_reviewers.docx`
   - `figure_thumbnails.pdf` (one page of small previews)
   - Cover note as email body or `v2_change_summary.md`.

4. **Send** (user task — handles the actual email).

5. **Co-author feedback window** ~1-2 weeks.

6. **Integration round.** Consolidate co-author edits; address comments;
   final pass. Target small enough for no second circulation.

7. **Pre-submission final-mile.**
   - All `[AUTHOR INPUT NEEDED]` markers resolved.
   - Reviewer affiliations and emails confirmed.
   - Diamond Apr 2026 proposal number in acknowledgments.
   - Data availability statement updated.
   - Final read-through.

**Deliverable.** v2 .docx bundle sent to co-authors → revised
manuscript ready for JAC submission.

## 10. Outstanding `[AUTHOR INPUT NEEDED]` items

1. **Beamtime wavelength** (Phase 1 step 1). Confirm λ for the
   24 Apr 2026 I11 MAC configuration.
2. **NWA 7465 metadata.** Lit shock stage, modal mineralogy, primary
   reference.
3. **NWA 5478 metadata.** Lit shock stage (S5? S6?), modal mineralogy,
   reference.
4. **NWA 7831 metadata.** Lit shock stage, opx/olv ratio, reference.
5. **NWA 6693 metadata.** Confirm classification (ungrouped achondrite,
   CC-paired with NWA 6704), shock stage (low — S1 or S2?), primary
   reference (Warren et al.?).
6. **Diamond proposal number** for the Apr 2026 beamtime (for
   acknowledgments and data availability).
7. **Decision on Fig 7 vs Fig S2.** Whether the capillary
   before/after demo goes in main text or supplementary, and
   which sample to demonstrate it with.
8. **Capillary subtraction citation.** Whether a reference is
   warranted (Billinge group, Wright et al.) and which one.
9. **Reviewer emails and affiliations** (still outstanding from
   the previous round): Ungár, Billinge, Langenhorst, Tsuchiyama,
   Toby + 5 alternates.

## 11. Risks and abort conditions

| Risk | Phase | Abort condition |
|---|---|---|
| Wavelength unknown | 1.1 | Halt until proposal/log book consulted. |
| Capillary subtraction disrupts PDF FWHM trend | 1.7 | Fall back to Option A (capillary as supplementary methods only) or escalate to user. |
| NWA 5478 Plag/Pyx ratio unexpectedly > 1 | 2.8 | Investigate before generating figures; may require petrographic re-examination. |
| Talampaya synchrotron W-A column length unexpectedly smaller than its Lab Co value (310 Å) | 2.8 | Investigate processing or measurement; may indicate sample-side issue. |
| Existing PDF first-peak FWHM trend with shock is broken by capillary subtraction | 1.6 / 2.4 | Reinterpret as "trend was partly capillary contribution" — narrative-disrupting; halt and re-scope manuscript story. |
| Co-author pushback on the speed of v2 round | 5 | Extend feedback window; offer second round if needed. |
| Submission timing pressure (external) | 5 | Negotiate; or fall back to "submit N=29 now, revise in review" with a quick-pivot cover letter. |

## 12. Success criteria

- v2 manuscript correctly cites N=33 and 8 parent-body groups throughout.
- §4.3 no longer contains a "would enable" sentence about the capillary.
- §3.3 Plag/Pyx reversal is anchored by ≥3 high-shock eucrites.
- All data-bearing figures (Figs 1-7) regenerated; Inkscape rework
  punch list issued; `submission/figures/` not modified without
  explicit author approval.
- PDF FWHM monotonic-with-shock trend (§3.x, Fig 3 or 6) survives
  capillary subtraction.
- Co-author cover note one page; v2 .docx track-changes intact.
- All `[AUTHOR INPUT NEEDED]` markers tracked, with v2 having
  placeholders acceptable and final submission having all resolved.

## 13. Files modified (summary)

### 13.1 Created

- `Llunr/n33_expansion/` — new analysis directory:
  - `phase1_validation/` — diagnostic figure + decision note
  - `data/i11_2026Apr_reduced/` — reduced .xy from .dat
  - `data/i11_corrected/` — capillary-subtracted .xy (new and existing)
  - `results/per_sample/` — per-sample W-A and PDF results (6 new analyses)
  - `results/per_sample_pdf_refits/` — re-fit PDFs for 18 existing I11 samples
  - `results/cross_instrument_pairs.csv` — Lab Co + I11 for NWA 5751, Talampaya
  - `scripts/` — analysis driver scripts (Python)
- `Paper1_JAC/survey_results_33samples.csv`
- `Paper1_JAC/paper1_figures/` (regenerated PNGs, several files)
- `Paper1_JAC/manuscript_paper1_JAC_v2_trackchanges.docx`
- `Paper1_JAC/manuscript_paper1_JAC_v2_clean.docx`
- `Paper1_JAC/submission/cover_letter_v2.docx`
- `Paper1_JAC/submission/v2_change_summary.md`

### 13.2 Modified

- `Paper1_JAC/manuscript_paper1_JAC.md`
- `Paper1_JAC/table1_sample_inventory.md`
- `Paper1_JAC/submission/cover_letter.md`
- `Paper1_JAC/submission/suggested_reviewers.md` (light)
- `Paper1_JAC/submission/build_auxiliary_docx.py` (cover letter body)

### 13.3 Untouched

- `xrd_profile/` package (any version): API consumer only.
- `Paper1_JAC/submission/figures/` (Inkscape hand-edited).
- `Paper1_JAC/survey_results_29samples.csv` (kept as the v1 reproducibility anchor).
- `Paper1_JAC/manuscript_paper1_JAC_submission_EAC may 9 2026-PJAM-MI.docx`
  (kept as the v1 approved snapshot).
- Existing supplementary figures and Table S1 from the Kleppe-response
  work (2026-05-11).
