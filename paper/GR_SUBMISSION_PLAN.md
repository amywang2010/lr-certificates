# Genome Research Spatial Omics SI — Compliance Map and Execution Plan (2026-09-08)

Primary venue locked: GR "Spatial Omics" special issue (deadline 2026-11-01).
Pivot on desk rejection: RECOMB 2027 (historical full-paper deadline ~Nov 20;
CFP pending). Working submission target: **Sep 29, 2026**, so the GR desk screen
(1-4 weeks) resolves before the RECOMB window closes. If GR passes desk, we are
under review and RECOMB is foregone by rule (no dual submission); that is the
success path.

## 1. What GR wants (sources: SI call text, read 2026-09-08; CSHL author-guideline
summaries; recent GR computational papers: spRefine gr.281001.125, CCC dissection
35:1400-1414, Nullstrap-style calibration work at partner venues)

| GR requirement / preference | Our status | Action |
|---|---|---|
| Scope: computational/statistical methods for spatial data analysis; cell-cell communication; cancer | direct textual match; editors Leslie/Ma/Moffitt are computational | keep; cite the SI themes in cover letter |
| Biology-first center of gravity. Their most-cited desk-reject pattern: "biology as demonstration dataset for a method"; methods-first abstract + benchmark-style figures flagged | **our biggest risk.** Current abstract opens "We introduce a certificate framework" | Reframe: abstract and Results lead with the certified biological findings; framework is the instrument |
| Abstract <= 250 words | 318 (violates) | rewrite to <= 250, biology-first |
| Running title <= 50 chars | 78 (violates) | shorten |
| Main text ~10,000-word cap | ~6,300 words | compliant; room for addition results |
| Figures <= 6, tables <= 2 in main text | 4 figures; verify all cited in text (only 2 distinct "Figure N" refs found) | audit citations; move any excess tables to supplement |
| Code + data availability mandatory, reviewer-ready at submission | repo public (github amywang2010/b03-certificates), hashed manifests, public 10x data with URLs recorded in the pre-registrations | reformat into GR-style Data Availability + Code Availability statements |
| CRediT author contributions, ORCID, COI, funding, 3-5 suggested reviewers (outside author institutions) | not yet prepared | submission package task |
| Reproducibility culture (referees "notably rigorous on code availability and benchmark comparisons") | strong: pre-registrations, manifests, referee re-derivations | present without advertising process |

## 2. The reframing plan (the decisive edit)

Title stays methodological (it is accurate), but the abstract and Figure 1 lead
with findings:
- Opening move: which ligand-receptor claims in tumor tissue survive segmentation
  uncertainty; currently most cannot be certified at all.
- Findings first: 19/60 breast contrasts selected as fragile collapse to
  non-identifiable (including pairs significant by conventional tests); ESR1-PGR
  mutual exclusion certified in all three lung regions and fully replicating under
  an independent segmenter; CDH1-EGFR autocrine certificate in lung mirroring
  CDH1-ERBB2 in breast; the cross-tissue fragility contrast on the same statistic.
- Framework second: exact max-flow bounds over a declared uncertainty set, matched
  null, pre-registered gates, minutes of compute.
- Results reordering: Section 5 (biological certificates) before the calibration
  results; framework details stay in Methods-style sections.

## 3. Limitations referees will demand fixed -> immediate fixes

| Limitation | Fix | Cost | Status |
|---|---|---|---|
| L5: B=1000 p-floor | rerun headline rows at B=10,000 (both tissues) | minutes-hours; one-page prereg addendum first | scheduled |
| L3: co-expression only | neighbor-effect statistic certificate demo (inside the frozen topic scope) | ~1 day | scheduled |
| L4: one independent segmenter | coverage + sign check under Xenium multimodal cell segmentation (public 10x datasets: lymph node 5K, lung preview) | ~2-3 days incl. prereg | scheduled |
| L1: single section per tissue | screen BioStudies S-BIAD2146 (31 human FFPE Xenium sections, 16 tissues) against the frozen pair list; add ONE section, vendor battery + matched null (proseg optional) | ~5-7 days; KILLABLE if panel overlap fails -> document, drop, scope L1 as remaining | scheduled with kill switch |
| L2/L6/L7: label freezing, kernel, crop | no changes; already scoped with sensitivities; L6 upside framing retained | none | hold |

Discipline: every fix pre-registered in an addendum before compute
numbering if deviations occur); same gates; no softening; a tissue that certifies
nothing is reported as such.

## 4. Dated execution plan (target: submit Sep 29)

- **Sep 8-9:** manuscript reframe (abstract <=250, running title, Results reorder,
  figure-citation audit, GR statements); start S-BIAD2146 screening and dataset
  downloads in parallel (network time hides behind writing time).
- **Sep 10-11:** B=10k rerun + neighbor-effect prereg and implementation.
- **Sep 12-16:** third-segmenter leg (prereg, run, report); fourth-tissue pipeline
  if screening passed.
- **Sep 17-19:** integration: new Results subsections, limitations update,
  supplement, figures (<=6).
- **Sep 20-22:** strict-reviewer audit pass v2; style sweep; number-by-number
  verification against new artifacts.
- **Sep 23-24:** bioRxiv deposit (GR permits simultaneous deposit).
- **Sep 25-29:** GR package: cover letter (SI keyword, editors, biology-first
  significance), CRediT, ORCID, COI, funding, data/code statements, 3-5 suggested
  reviewers, format conversion, submit.
- **Buffer:** Sep 29-Oct 14 spare; even a Sep 29 slip keeps desk decision before
  mid-November.

## 5. Cover-letter skeleton (to draft Sep 25)

- SI keyword "Spatial Omics Special Issue"; addressed to Leslie/Ma/Moffitt.
- Paragraph 1: the reader's question (does this conclusion survive a different
  defensible segmentation?) and the biological answer we deliver.
- Paragraph 2: findings (fragile-claim collapse; certified exclusion atlas;
  cross-tissue contrast) in two sentences, no overclaim.
- Paragraph 3: what the method is (exact bounds, matched null, minutes of
  compute) and why the SI scope fits, quoting one call phrase.
- Paragraph 4: confirmations: not under consideration elsewhere; bioRxiv DOI;
  code/data reviewer-ready; all authors approve.
