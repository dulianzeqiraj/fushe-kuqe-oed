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

Link the GitHub account at <https://zenodo.org/account/settings/github/>,
switch the repository on, then cut a release:

```bash
git tag -a v1.0.0 -m "Version accompanying the submitted manuscript"
git push origin v1.0.0
```

Zenodo mints a DOI from the release. `.zenodo.json` supplies the title,
licence, keywords and the four related identifiers.

## 4. Put the real URLs into the manuscript

The two links are the only values in the paper that do not come out of a run,
so they are supplied from the environment rather than edited into the text:

```bash
cd "C:/Users/d_zeq/OneDrive/Desktop/PROJEKTE 2 ARTIKUJ/MREKULLIA/Q1_RIPUNIM_2026/code"
FK_GITHUB_URL=https://github.com/<user>/fushe-kuqe-oed \
FK_ZENODO_DOI=https://doi.org/10.5281/zenodo.<id> \
python fill_manuscript.py && python build_docx.py && python verify_docx.py
```

`verify_docx.py` refuses the document if a placeholder survived, if anything
other than the MATLAB listings follows the references, if an em dash appears,
or if any listing has drifted from the file on disk.

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
