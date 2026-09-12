# What changed from the earlier draft, and why

The file this package replaces is

`MBAS REFERNCAVE DUHEN VETEM KODET MATLAB  FUSHE_KUQE_BAYESIAN_ARTICLE_IMPROVED_Q1.docx`

It has not been modified. It is used only as the style reference for the new
document, so the page setup, fonts and named paragraph styles carry over.

Its name asked for the MATLAB code to follow the references. It did, through
v1.0.2. From v1.0.3 the code is in the repository only: the manuscript is
twenty-one pages of paper rather than a hundred of listing, and the Data
Availability section says where the pipeline is. `fk_verify_docx.m` now fails
if a code listing reappears in the document.

## Why the draft could not be polished into a submission

Five things were found on inspection, each checked rather than assumed.

**1. The numerical results were synthetic.** `MATLAB_CODES_COMPLETE_PACKAGE/generate_synthetic_data.m`
states in its own header that it generates data "matching the statistical
properties reported in the manuscript", with the manuscript's numbers hard
coded (`pathway1_mean = 0.285`, `fused_CV = 9.9`) and `rng(42)` well positions.
`sec63_matlab.m` builds 37 wells from `randn` in UTM zone 34N. The real
workbook is in Gauss-Krüger zone 4 on the Krasovsky 1940 datum and holds 180
wells. Figures produced this way cannot support a field claim.

**2. Theorem 1 and its numbers are already published.** The Dual-Pathway
fusion, the Cramér-Rao optimality proof, `e_F = 1.00`, the 33 per cent
uncertainty reduction, `I₁ = 566.89`, `I₂ = 384.47`, `I_total = 951.36` and
`n_e = 0.282 ± 0.032` all appear in Zeqiraj et al. (2026),
*Journal of Hazardous Materials Advances* 23, 101261
(doi:10.1016/j.hazadv.2026.101261). Resubmitting them as new results is
duplicate publication.

**3. The independence assumption is contradicted by your own later paper.**
Zeqiraj and Beqiraj (2026), *Journal of Contaminant Hydrology* 283, 105086
(doi:10.1016/j.jconhyd.2026.105086) is titled "when the errors are not
independent" and shows that the Vukovic-Soro and Kozeny-Carman pathways share
grain-size information, with facies correlations of 0.15, 0.35 and 0.58.
Section 3.6 of the draft verifies independence on the same two pathways.

**4. The abstract reused a published sentence almost verbatim.** The draft's
closing sentence and the closing sentence of the JHMA abstract differ in a few
words. Crossref Similarity Check would flag it.

**5. Table 4 was empty.** Every cell read `[X]`, the authors were
`[Author 1]` to `[Author 3]`, and the repository was "to be provided upon
acceptance".

## What the new paper is

The parts of the draft that were genuinely unpublished are the information
geometry, the implicit particle machinery and the monitoring-network design.
Those became the paper. The fusion is now an input, cited to the two published
papers and not re-derived, which also removes the contradiction in point 3
because the correlation-aware variance is the one used.

The contribution is now: what the existing network is worth, measured in
degrees of freedom for signal; why a steady concentration field cannot
constrain effective porosity at all; a two-line proof that weighting the
parameter space leaves D-optimality unchanged; a goal-oriented criterion that
does change the design; and a computed answer to the concentration-against-
pumping-test question, which comes out inverted.

## What is different in the data

| | draft | this paper |
|---|---|---|
| source | synthetic, `rng(42)` | `Fushe_Kuqe_All_Data.xlsx` (Desktop copy) |
| wells | 37 | 180 total, 43 with K, 31 with NO₃ |
| area | 1,200 km² | 349 km², the digitised polygon |
| coordinates | UTM 34N | Gauss-Krüger zone 4, Krasovsky 1940 |
| flow direction | SE→NW rotating to E→W | NE→SW, from Cenameri and Beqiraj (2016) |
| second site | none | 44 USGS wells, Pugh (2023) |

The flow direction changed because the three seawater-intrusion sheets in the
Desktop copy of the workbook carry the primary description
(doi:10.12681/bgsg.11772), which gives NE→SW with the wedge advancing the
other way. The earlier rotating field does not come from that source.

## Two things in the workbook worth knowing

**The lithological label does not separate K.** One-way ANOVA of hydraulic
conductivity on the coarse/medium/fine gravel label over the 43 wells gives
p = 0.44. Fine gravel includes three wells at 200 m d⁻¹; coarse gravel
includes one at 50 m d⁻¹. Facies-keyed parameter assignment is therefore not
supported by this file, and the prior in the new paper is keyed on measured
log₁₀K instead.

**The nitrate field is not spatially predictable at this density.** Leave-one-
out cross-validation over the 24 nitrate wells inside the aquifer: ordinary
kriging R² between −0.11 and −0.28 at every correlation length from 1 to
16 km, inverse-distance weighting down to −0.74, DRASTIC alone −0.19. The only
thing that beats the network mean is a linear trend in easting and northing, at
R² = 0.16. This is reported as a result, and it is why the paper's design
argument is about what to measure next rather than about fitting what exists.

## Verification built into the package

- `run_all.m` writes `results/results.json`; nothing else may be quoted.
- The whole chain is MATLAB. The document is written as WordprocessingML
  by `fk_md2docx.m`, with the equations converted to Office MathML by
  `fk_latex2omml.m`, so they stay editable in Word. There is no Python and
  no pandoc.
- `fk_fill_manuscript.m` refuses to produce the manuscript if any `{{placeholder}}`
  has no value in that file, and lists values computed but never cited.
- `fk_verify_docx.m` checks that nothing but the MATLAB listings follows the
  references, that no placeholder survived, that there is no em dash, and that
  each listing's SHA-256 prefix matches the file on disk.
- `fk_sensitivity.m` compares its analytic Jacobian against finite differences
  at three step sizes and reports the convergence order.
- `mrva_replicate.m` recomputes every DRASTIC index from its seven ratings and
  stops the run if any fails to reproduce.
- Every reference DOI was resolved against Crossref before it entered the list.
  Two entries have no DOI: Fedorov (1972), verified in Open Library, and
  Krause et al. (2008), verified at jmlr.org.

## One defect that reached a delivered document, and the check that now catches it

The first build of this manuscript, committed as eecb296 and pushed, contained
`$lpha_L$` where it should have read an equation. A stray byte had entered the
markdown during an edit: `lpha` became a BEL control character followed by
"lpha". The converter in use at the time dropped the byte without complaining
and the text reached print as a word.

Two things changed because of it. `fk_fill_manuscript.m` now refuses to write a
manuscript containing any control character, and `fk_verify_docx.m` checks the
finished document for control characters and for stray dollar signs, the latter
being what an unconverted equation leaves behind. Both checks are in the
listings; neither existed when the defect shipped.

## The figures were missing from the document, and now are not

The first builds carried the figure captions and none of the figures. The
document said nothing about it: captions read normally, and the count of
tables and equations was right, so every check passed. He caught it.

The document writer now embeds pictures. `fk_md2docx.m` understands
`![alt](path)` and writes the drawing markup, adds the image parts, amends the
relationships and declares the PNG content type; `fk_png_size.m` reads the
dimensions so a figure is set to the text measure with its aspect preserved;
`fk_zip_write.m` replaces the old single-part writer so that parts can be added
as well as replaced. `fk_verify_docx.m` now counts pictures against captions
and fails if any caption has none.

The figures sit in a Figures section after the CRediT statement and before the
references, because after the references the document carries nothing but the
MATLAB listings.

## Provenance of manuscript/template.docx

The template carries the styles, numbering, theme and page setup and no content.
It was made by `fk_make_template.m` from the first build of this manuscript,
which in turn took its styles from the author's original draft. Nothing in the
scientific content passes through it.

## Two dates to be aware of

Crossref registers Cenameri and Beqiraj as 2017 although the volume is the
2016 congress proceedings and the workbook cites 2016; the reference list uses
2016. Crossref registers Pugh as 2025 although the report is numbered
2023-5101; the reference list uses 2023 with the report number.
