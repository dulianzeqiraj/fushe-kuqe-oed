function fk_verify_docx(docx)
%FK_VERIFY_DOCX  Check the submission document before it goes anywhere.
%
%  Six checks, each of which has caught something at least once:
%
%    1. Nothing but the MATLAB listings follows the reference list.  That is
%       the one structural requirement the document has.
%    2. No unresolved {{placeholder}} survived into the file.
%    3. No em dash anywhere.
%    4. No control character anywhere.  This check exists because one reached
%       a delivered document on 12 September 2026: a stray byte turned
%       $\alpha_L$ into the word "lpha", and the converter dropped it in
%       silence.
%    5. No stray dollar sign, which is what an unconverted equation leaves.
%    6. Every listing in the appendix matches the file on disk, compared by
%       the SHA-256 prefix printed beside it.
%    7. Every figure caption has a picture above it.  This check exists
%       because the first builds carried the captions and none of the
%       figures, and nothing in the document said so.

if nargin < 1
    cfg = fk_config();
    docx = fullfile(cfg.dir_code, '..', 'manuscript', ...
                    'Zeqiraj_FusheKuqe_OED_manuscript.docx');
end
if ~exist(docx, 'file')
    error('fk_verify_docx:missing', 'no document at %s', docx);
end

[paras, styles, text] = local_paragraphs(docx);
fails = {};

% ---- 1. structure ---------------------------------------------------------
iref = find(strcmpi(strtrim(paras), 'references'), 1);
iapp = find(startsWith(lower(strtrim(paras)), 'appendix a'), 1);
if isempty(iref)
    fails{end+1} = 'no References heading';
elseif isempty(iapp)
    fails{end+1} = 'no Appendix A heading';
elseif iapp < iref
    fails{end+1} = 'the appendix precedes the references';
else
    fprintf(['[verify] %d paragraphs, references at %d, appendix at %d, ' ...
             '%d after it\n'], numel(paras), iref, iapp, ...
            numel(paras) - iapp);
end

% The prose checks below apply to the paper, not to the code listings, which
% legitimately contain dollar signs and any character a MATLAB file may hold.
if ~isempty(iapp)
    prose = strjoin(paras(1:iapp), newline);
else
    prose = text;
end

% ---- 2. placeholders ------------------------------------------------------
ph = regexp(prose, '\{\{\w+\}\}', 'match');
if ~isempty(ph)
    fails{end+1} = ['unresolved placeholders: ' strjoin(unique(ph), ', ')];
end

% ---- 3. em dashes ---------------------------------------------------------
nem = sum(double(prose) == 8212);
if nem > 0
    fails{end+1} = sprintf('%d em dash(es)', nem);
end

% ---- 4. control characters ------------------------------------------------
ctrl = double(text) < 32 & ~ismember(double(text), [9 10 13]);
if any(ctrl)
    j = find(ctrl, 1);
    fails{end+1} = sprintf('control character U+%04X near "%s"', ...
        double(text(j)), text(max(1, j-40):min(numel(text), j+20)));
end

% ---- 5. stray dollar signs ------------------------------------------------
% Equations become OMML, which carries no dollar sign; one left in the text
% means an equation was not converted.
nd = sum(double(prose) == 36);
if nd > 0
    j = find(double(prose) == 36, 1);
    fails{end+1} = sprintf('%d dollar sign(s), first near "%s"', nd, ...
        prose(max(1, j-50):min(numel(prose), j+30)));
end

% ---- 6. listings match the code -------------------------------------------
cfg = fk_config();
tok = regexp(text, '([A-Za-z0-9_]+\.m)\s+(\d+)\s+([0-9a-f]{16})', 'tokens');
nchk = 0;
for k = 1:numel(tok)
    name = tok{k}{1};
    dg = tok{k}{3};
    path = fullfile(cfg.dir_code, name);
    if ~exist(path, 'file')
        fails{end+1} = ['listed file not on disk: ' name]; %#ok<AGROW>
        continue
    end
    real = fk_sha256(path);
    if ~strcmp(real, dg)
        fails{end+1} = sprintf('listing out of date: %s (document %s, disk %s)', ...
                               name, dg, real); %#ok<AGROW>
    end
    nchk = nchk + 1;
end
fprintf('[verify] %d listing checksums checked\n', nchk);
if nchk == 0
    fails{end+1} = 'no listing checksums found';
end

% ---- 7. figures present ---------------------------------------------------
raw = local_rawxml(docx);
nimg = numel(strfind(raw, '<w:drawing>'));
ncap = sum(cellfun(@(t) ~isempty(regexp(strtrim(t), ...
           '^(Figure \d+\.|Graphical abstract\.)', 'once')), paras));
fprintf('[verify] %d figures embedded, %d captions\n', nimg, ncap);
if ncap == 0
    fails{end+1} = 'no figure captions found';
elseif nimg < ncap
    fails{end+1} = sprintf('%d captions but only %d figures embedded', ...
                           ncap, nimg);
end

% ---- report ---------------------------------------------------------------
if ~isempty(fails)
    fprintf('\n[verify] FAILED\n');
    for k = 1:numel(fails)
        fprintf('  - %s\n', fails{k});
    end
    error('fk_verify_docx:failed', '%d check(s) failed', numel(fails));
end
fprintf('[verify] all checks passed\n');
end

% =========================================================================
function [paras, styles, text] = local_paragraphs(docx)
tmp = tempname;
mkdir(tmp);
c = onCleanup(@() rmdir(tmp, 's')); %#ok<NASGU>
unzip(docx, tmp);
x = fileread(fullfile(tmp, 'word', 'document.xml'));

pieces = regexp(x, '<w:p[ >].*?</w:p>|<w:p/>', 'match');
paras = cell(numel(pieces), 1);
styles = cell(numel(pieces), 1);
for k = 1:numel(pieces)
    s = pieces{k};
    t = regexp(s, '<(?:w|m):t[^>]*>(.*?)</(?:w|m):t>', 'tokens');
    txt = '';
    for j = 1:numel(t)
        txt = [txt local_unesc(t{j}{1})]; %#ok<AGROW>
    end
    paras{k} = txt;
    st = regexp(s, '<w:pStyle w:val="([^"]+)"', 'tokens', 'once');
    if isempty(st), styles{k} = ''; else, styles{k} = st{1}; end
end
text = strjoin(paras, newline);
end

function x = local_rawxml(docx)
tmp = tempname;
mkdir(tmp);
c = onCleanup(@() rmdir(tmp, 's')); %#ok<NASGU>
unzip(docx, tmp);
x = fileread(fullfile(tmp, 'word', 'document.xml'));
end

function s = local_unesc(s)
s = strrep(s, '&lt;', '<');
s = strrep(s, '&gt;', '>');
s = strrep(s, '&quot;', '"');
s = strrep(s, '&apos;', '''');
s = strrep(s, '&amp;', '&');
end
