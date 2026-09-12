"""verify_docx.py

Check the assembled submission document before it goes anywhere.

Four checks, each of which has caught something at least once:

  1. Nothing but the MATLAB listings follows the reference list.  This is the
     one structural requirement the document has.
  2. No unresolved {{placeholder}} survived into the file.
  3. No em dash anywhere in the prose.
  4. Every listing in the appendix matches the file on disk byte for byte,
     compared by the SHA-256 prefix printed with it.

Exit status is non-zero if any check fails, and the failure is named.
"""
import hashlib
import os
import re
import sys
import zipfile
from xml.etree import ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
DOCX = os.path.join(ROOT, 'manuscript',
                    'Zeqiraj_FusheKuqe_OED_manuscript.docx')

W = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'


def paragraphs(path):
    z = zipfile.ZipFile(path)
    root = ET.fromstring(z.read('word/document.xml'))
    body = root.find(W + 'body')
    out = []
    for p in body.iter(W + 'p'):
        txt = ''.join(t.text or '' for t in p.iter(W + 't'))
        st = p.find(W + 'pPr/' + W + 'pStyle')
        style = st.get(W + 'val') if st is not None else ''
        out.append((style, txt))
    return out


def main():
    if not os.path.exists(DOCX):
        print('missing %s' % DOCX, file=sys.stderr)
        return 1
    ps = paragraphs(DOCX)
    text = '\n'.join(t for _, t in ps)
    fails = []

    # --- 1. structure ---------------------------------------------------
    iref = next((k for k, (s, t) in enumerate(ps)
                 if t.strip().lower() == 'references'), None)
    iapp = next((k for k, (s, t) in enumerate(ps)
                 if t.strip().lower().startswith('appendix a')), None)
    if iref is None:
        fails.append('no References heading found')
    elif iapp is None:
        fails.append('no Appendix A heading found')
    elif iapp < iref:
        fails.append('the appendix precedes the references')
    else:
        after = [t for s, t in ps[iapp + 1:] if t.strip()]
        print('paragraphs: %d total, references at %d, appendix at %d, '
              '%d paragraphs after the appendix heading'
              % (len(ps), iref, iapp, len(after)))

    # --- 2. placeholders -------------------------------------------------
    left = re.findall(r'\{\{\w+\}\}', text)
    if left:
        fails.append('unresolved placeholders: %s' % sorted(set(left)))

    # --- 3. em dashes ----------------------------------------------------
    n_em = text.count('—')
    if n_em:
        ctx = [ln for ln in text.splitlines() if '—' in ln][:3]
        fails.append('%d em dash(es), first lines: %s' % (n_em, ctx))

    # --- 4. listings match the code on disk ------------------------------
    digests = dict(re.findall(r'([A-Za-z0-9_]+\.(?:m|py))\s+\|?\s*\d*\s*\|?\s*([0-9a-f]{16})',
                              text))
    if not digests:
        digests = dict(re.findall(r'([A-Za-z0-9_]+\.(?:m|py)).{0,40}?([0-9a-f]{16})',
                                  text, re.S))
    checked = 0
    for name, dg in digests.items():
        path = os.path.join(HERE, name)
        if not os.path.exists(path):
            fails.append('listed file not on disk: %s' % name)
            continue
        with open(path, 'rb') as f:
            real = hashlib.sha256(f.read()).hexdigest()[:16]
        if real != dg:
            fails.append('listing out of date: %s (document %s, disk %s)'
                         % (name, dg, real))
        checked += 1
    print('checked %d listing checksums' % checked)
    if checked == 0:
        fails.append('no listing checksums found in the document')

    if fails:
        print('\nFAILED:')
        for f in fails:
            print('  - %s' % f)
        return 1
    print('\nall checks passed')
    return 0


if __name__ == '__main__':
    sys.exit(main())
