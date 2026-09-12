# Publishing the package

Nothing here has been pushed anywhere. These are the steps, in order, and the
two decisions that have to be made before the first one.

## Before anything is pushed

**Do you have the right to release the well data?** The workbook is a
compilation from Albanian Geological Survey and IGJEUM records. The repository
puts the coordinates, hydraulic conductivities, nitrate concentrations and
DRASTIC indices of 180 monitoring points into the public domain permanently.
If either body has to agree, it has to agree first; a Zenodo deposit cannot be
withdrawn once it has a DOI. If the answer is no, the alternative is to release
the code with the derived quantities only, and to state in the Data
Availability section that the raw well table is available from the
corresponding author on request. The manuscript currently promises the full
release.

**ORCID.** Both files now carry 0000-0002-5305-3048. That was not supplied; it was found in the
Crossref record for doi:10.1016/j.rineng.2021.100242, on dulianzeqiraj.github.io, and confirmed
against the ORCID public API, which returns the name Dulian Zeqiraj. Check it before deposit.

## 1. Local repository

```bash
cd "C:/Users/d_zeq/OneDrive/Desktop/PROJEKTE 2 ARTIKUJ/MREKULLIA/Q1_RIPUNIM_2026"
git init -b main
git add .
git commit -m "Information-geometric monitoring network design for the Fushe-Kuqe aquifer"
```

`.gitignore` already excludes `results/run_all.mat` and `results/stage_*.mat`,
which are large and regenerable. `results/results.json` is committed, because
the manuscript quotes from it and a reader has to be able to check.

## 2. GitHub

```bash
gh repo create fushe-kuqe-oed --public --source=. --remote=origin --push
```

or, without the CLI, create the repository on github.com and then

```bash
git remote add origin https://github.com/dulianzeqiraj/fushe-kuqe-oed.git
git push -u origin main
```

## 3. Zenodo

Done on 12 September 2026. The archive is

- concept DOI, always the latest version: <https://doi.org/10.5281/zenodo.22726726>
- version DOI for v1.0.0: <https://doi.org/10.5281/zenodo.22726727>

The manuscript cites the concept DOI, so it stays valid as versions are added.

One thing to know if this is ever repeated on a new repository. Zenodo archives
only releases made after the repository is switched on, and its repository list
comes from a cached sync, so a repository created the same day will not appear
until Sync now is pressed. The first three releases here produced nothing
because the Zenodo account was still linked to an older GitHub account, and it
could not see this repository at all.

To cut a further version: link the GitHub account at
<https://zenodo.org/account/settings/github/>, check the repository is on, then

```bash
git tag -a v1.0.0 -m "Version accompanying the submitted manuscript"
git push origin v1.0.0
```

Zenodo mints a DOI from the release. `.zenodo.json` supplies the title,
licence, keywords and the four related identifiers.

## 4. Put the real URLs into the manuscript

The two links are the only values in the paper that do not come out of a run,
so they are passed as arguments rather than edited into the text. The GitHub
URL is already the default; only the Zenodo DOI has to be supplied:

```bash
cd "C:/Users/d_zeq/OneDrive/Desktop/PROJEKTE 2 ARTIKUJ/MREKULLIA/Q1_RIPUNIM_2026/code" && matlab -batch "fk_fill_manuscript('', 'https://doi.org/10.5281/zenodo.NRI'); fk_build_docx; fk_verify_docx"
```

`fk_verify_docx` refuses the document if a placeholder survived, if anything
other than the MATLAB listings follows the references, if an em dash, a control
character or a stray dollar sign appears, or if any listing has drifted from
the file on disk.

## 5. What to submit

| File | What it is |
|---|---|
| `manuscript/Zeqiraj_FusheKuqe_OED_manuscript.docx` | the manuscript, listings after the references |
| `manuscript/Highlights.docx` | the five highlights |
| `manuscript/Cover_Letter.docx` | the cover letter |
| `figures/fig1..fig7.tif` | figures at 600 dpi |
| `figures/graphical_abstract.tif` | graphical abstract |

## One practical note

The MATLAB licence on this machine reported five days to expiry on 12 September
2026. The pipeline needs it to reproduce anything.
