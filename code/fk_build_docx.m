function fk_build_docx()
%FK_BUILD_DOCX  Assemble the submission documents.
%
%  Produces, in manuscript/:
%    Zeqiraj_FusheKuqe_OED_manuscript.docx   the paper, MATLAB listings after
%                                            the references and nothing else
%    Highlights.docx
%    Cover_Letter.docx
%
%  The listings are read from the .m files themselves, so the appendix cannot
%  drift from the code that produced the numbers, and each is identified by a
%  SHA-256 prefix that fk_verify_docx checks against the file on disk.
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

% ------------------------------------------------- the listings, in order --
order = {'fk_config.m', 'fk_load_data.m', 'fk_predictability.m', ...
         'fk_build_grid.m', 'fk_prior_field.m', 'fk_forward.m', ...
         'fk_calibrate.m', 'fk_sensitivity.m', 'fk_design.m', ...
         'oed_local.m', 'fk_particles.m', 'mrva_replicate.m', ...
         'fk_sweep.m', 'fk_write_results.m', 'fk_figures.m', ...
         'fk_prepare_data.m', 'fk_fill_manuscript.m', 'fk_latex2omml.m', ...
         'fk_xmlesc.m', 'fk_check_omml.m', 'fk_md2docx.m', ...
         'fk_zip_replace.m', 'fk_make_template.m', 'fk_build_docx.m', ...
         'fk_verify_docx.m', 'fk_sha256.m', 'run_all.m'};

body = fileread(src);
parts = {strtrim(body), '', '', '## Appendix A. MATLAB code', '', ...
   ['The listings below are the complete pipeline, in the order it runs ' ...
    'them. They are reproduced verbatim from the released repository; the ' ...
    'sixteen-character SHA-256 prefix given with each file identifies the ' ...
    'exact version that produced every number in this paper. Nothing ' ...
    'follows the listings.'], ''};

tbl = {'| File | Lines | SHA-256 prefix |', '|---|---|---|'};
lst = {};
for k = 1:numel(order)
    path = fullfile(cfg.dir_code, order{k});
    if ~exist(path, 'file')
        error('fk_build_docx:listing', 'missing listing %s', order{k});
    end
    code = fileread(path);
    nl = numel(strfind(code, newline)) + 1;
    dg = fk_sha256(path);
    tbl{end+1} = sprintf('| `%s` | %d | `%s` |', order{k}, nl, dg); %#ok<AGROW>
    lst{end+1} = sprintf('### A.%d  `%s`', k, order{k}); %#ok<AGROW>
    lst{end+1} = ''; %#ok<AGROW>
    lst{end+1} = '```matlab'; %#ok<AGROW>
    lst{end+1} = strip_trailing(code); %#ok<AGROW>
    lst{end+1} = '```'; %#ok<AGROW>
    lst{end+1} = ''; %#ok<AGROW>
    fprintf('%-24s %5d lines  %s\n', order{k}, nl, dg);
end

full = strjoin([parts, tbl, {''}, lst], newline);
fullmd = fullfile(man, 'manuscript_with_code.md');
fid = fopen(fullmd, 'w', 'n', 'UTF-8');
fwrite(fid, unicode2native(full, 'UTF-8'));
fclose(fid);

fk_md2docx(full, tpl, out);

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
