"""build_docx.py

Assemble the submission document: the filled manuscript, then the references,
then nothing but the MATLAB listings.

The listings are read from the .m files themselves, so the appendix cannot
drift from the code that produced the numbers.  A checksum of each file is
printed, and the same checksums are written into the appendix header, so a
reader can verify that the listing matches the released repository.

Pandoc converts to .docx using the author's existing manuscript as the style
reference, which carries the page setup, fonts and named paragraph styles over
to the new document.  The original file is opened read-only and is never
modified.
"""
import hashlib
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
MAN = os.path.join(ROOT, 'manuscript')
SRC = os.path.join(MAN, 'manuscript_filled.md')
FULL = os.path.join(MAN, 'manuscript_with_code.md')
OUT = os.path.join(MAN, 'Zeqiraj_FusheKuqe_OED_manuscript.docx')

REFERENCE_DOCX = (r'C:/Users/d_zeq/OneDrive/Desktop/PROJEKTE 2 ARTIKUJ/MREKULLIA/'
                  r'MBAS REFERNCAVE DUHEN VETEM KODET MATLAB  '
                  r'FUSHE_KUQE_BAYESIAN_ARTICLE_IMPROVED_Q1.docx')

# The order the listings appear in, which is the order the pipeline runs them.
ORDER = [
    'prepare_data.py',
    'fk_config.m',
    'fk_load_data.m',
    'fk_predictability.m',
    'fk_build_grid.m',
    'fk_prior_field.m',
    'fk_forward.m',
    'fk_calibrate.m',
    'fk_sensitivity.m',
    'fk_design.m',
    'oed_local.m',
    'fk_particles.m',
    'mrva_replicate.m',
    'fk_sweep.m',
    'fk_write_results.m',
    'fk_figures.m',
    'run_all.m',
]


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        h.update(f.read())
    return h.hexdigest()[:16]


def main():
    if not os.path.exists(SRC):
        print('run fill_manuscript.py first: %s is missing' % SRC, file=sys.stderr)
        return 1

    with open(SRC, encoding='utf-8') as f:
        body = f.read()

    parts = [body.rstrip(), '', '', '## Appendix A. MATLAB code', '']
    parts.append(
        'The listings below are the complete pipeline, in the order it runs '
        'them. They are reproduced verbatim from the released repository; the '
        'sixteen-character SHA-256 prefix given with each file identifies the '
        'exact version that produced every number in this paper. Nothing '
        'follows the listings.')
    parts.append('')

    tbl = ['| File | Lines | SHA-256 prefix |', '|---|---|---|']
    listings = []
    for name in ORDER:
        path = os.path.join(HERE, name)
        if not os.path.exists(path):
            print('missing listing: %s' % name, file=sys.stderr)
            return 1
        with open(path, encoding='utf-8') as f:
            code = f.read()
        nl = code.count('\n') + 1
        dg = sha(path)
        tbl.append('| `%s` | %d | `%s` |' % (name, nl, dg))
        lang = 'python' if name.endswith('.py') else 'matlab'
        listings.append('### A.%d  `%s`' % (len(listings) + 1, name))
        listings.append('')
        listings.append('```%s' % lang)
        listings.append(code.rstrip())
        listings.append('```')
        listings.append('')
        print('%-24s %5d lines  %s' % (name, nl, dg))

    parts.extend(tbl)
    parts.append('')
    parts.extend(listings)

    with open(FULL, 'w', encoding='utf-8') as f:
        f.write('\n'.join(parts))
    print('\nwrote %s' % FULL)

    cmd = ['pandoc', FULL, '-f', 'markdown', '-t', 'docx',
           '--reference-doc', REFERENCE_DOCX,
           '--highlight-style', 'monochrome',
           '-o', OUT]
    print(' '.join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout)
        print(r.stderr, file=sys.stderr)
        return r.returncode
    print('wrote %s (%.0f kB)' % (OUT, os.path.getsize(OUT) / 1024))

    # --- the two ancillary files a submission needs, from the same sources ---
    hl = []
    take = False
    for line in body.splitlines():
        if line.startswith('## Highlights'):
            take = True
            hl.append('Highlights')
            hl.append('')
            continue
        if take and line.startswith('## '):
            break
        if take and line.strip():
            hl.append(line)
    hi_md = os.path.join(MAN, '_highlights.md')
    with open(hi_md, 'w', encoding='utf-8') as f:
        f.write('\n'.join(hl))
    for src, dst in ((hi_md, 'Highlights.docx'),
                     (os.path.join(MAN, 'cover_letter.md'), 'Cover_Letter.docx')):
        if not os.path.exists(src):
            continue
        c = ['pandoc', src, '-f', 'markdown', '-t', 'docx',
             '--reference-doc', REFERENCE_DOCX, '-o', os.path.join(MAN, dst)]
        rr = subprocess.run(c, capture_output=True, text=True)
        if rr.returncode == 0:
            print('wrote %s' % dst)
        else:
            print(rr.stderr, file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
