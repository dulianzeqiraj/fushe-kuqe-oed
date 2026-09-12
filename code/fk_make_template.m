function fk_make_template(src, dst)
%FK_MAKE_TEMPLATE  Strip a .docx down to its styles, to serve as a template.
%
%  fk_make_template(src, dst)
%
%  Keeps every part of src except the body of word/document.xml, which is
%  emptied apart from the section properties.  The result carries the styles,
%  numbering definitions, theme, fonts and page setup, and nothing else.
%
%  The template in manuscript/template.docx was made this way and is committed
%  so that a reader does not have to reproduce it; this function is here so
%  that they can.

if nargin < 2
    cfg = fk_config();
    man = fullfile(cfg.dir_code, '..', 'manuscript');
    if nargin < 1
        src = fullfile(man, 'Zeqiraj_FusheKuqe_OED_manuscript.docx');
    end
    dst = fullfile(man, 'template.docx');
end

tmp = tempname;
mkdir(tmp);
c = onCleanup(@() rmdir(tmp, 's')); %#ok<NASGU>
unzip(src, tmp);
d = fileread(fullfile(tmp, 'word', 'document.xml'));

hdr = regexp(d, '^.*?<w:body>', 'match', 'once');
if isempty(hdr)
    error('fk_make_template:body', '%s has no w:body', src);
end
sp = regexp(d, '<w:sectPr[^>]*>.*?</w:sectPr>', 'match', 'once');
if isempty(sp)
    sp = regexp(d, '<w:sectPr[^>]*/>', 'match', 'once');
end

xml = [hdr sp '</w:body></w:document>'];
fk_zip_replace(src, dst, 'word/document.xml', xml);

z = java.util.zip.ZipFile(dst);
np = 0;
en = z.entries();
while en.hasMoreElements(), en.nextElement(); np = np + 1; end
z.close();
fprintf('[template] %s written, %d parts, body emptied\n', dst, np);
end
