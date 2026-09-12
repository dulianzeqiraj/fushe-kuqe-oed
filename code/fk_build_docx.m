function fk_build_docx()
%FK_BUILD_DOCX  Assemble the submission documents.
%
%  Produces, in manuscript/:
%    Zeqiraj_FusheKuqe_OED_manuscript.docx   the paper, ending at the
%                                            reference list
%    Highlights.docx
%    Cover_Letter.docx
%
%  The MATLAB listings are not printed in the manuscript.  They live in the
%  repository cited in the Data Availability section, which is where a reader
%  can run them rather than only read them; fk_verify_docx checks that none
%  has crept back into the document.
%
%  Run fk_fill_manuscript first.

cfg = fk_config();
man = fullfile(cfg.dir_code, '..', 'manuscript');
src = fullfile(man, 'manuscript_filled.md');
tpl = fullfile(man, 'template.docx');
out = fullfile(man, 'Zeqiraj_FusheKuqe_OED_manuscript.docx');

if ~exist(src, 'file')
    error('fk_build_docx:missing', 'run fk_fill_manuscript first: %s', src);
end
if ~exist(tpl, 'file')
    error('fk_build_docx:template', 'no template: %s', tpl);
end

fk_check_omml(src);

body = fileread(src);

fk_md2docx(body, tpl, out);

% --------------------------------------------- highlights and cover letter --
hi = regexp(body, '## Highlights\s*(.*?)\n## ', 'tokens', 'once');
if ~isempty(hi)
    hmd = ['## Highlights' newline newline strtrim(hi{1}) newline];
    fk_md2docx(hmd, tpl, fullfile(man, 'Highlights.docx'));
end
cl = fullfile(man, 'cover_letter.md');
if exist(cl, 'file')
    fk_md2docx(fileread(cl), tpl, fullfile(man, 'Cover_Letter.docx'));
end
end

% =========================================================================
function s = strip_trailing(s)
while ~isempty(s) && (s(end) == newline || s(end) == char(13))
    s(end) = [];
end
end
