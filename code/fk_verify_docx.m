function fk_verify_docx(docx)
%FK_VERIFY_DOCX  Check the submission document before it goes anywhere.
%
%  Six checks, each of which has caught something at least once:
%
%    1. The reference list is the last thing in the document.  The MATLAB
%       code lives in the repository, not in the paper.
%    2. No unresolved {{placeholder}} survived into the file.
%    3. No em dash anywhere.
%    4. No control character anywhere.  This check exists because one reached
%       a delivered document on 12 September 2026: a stray byte turned
%       $\alpha_L$ into the word "lpha", and the converter dropped it in
%       silence.
%    5. No stray dollar sign, which is what an unconverted equation leaves.
%    6. No code listing is present, so the pipeline cannot creep back into
%       the manuscript unnoticed.
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
% The reference entries follow the References heading, as they must.  What
% must not follow is another section: the MATLAB code used to be printed here
% and now lives in the repository instead.
iref = find(strcmpi(strtrim(paras), 'references'), 1);
if isempty(iref)
    fails{end+1} = 'no References heading';
else
    later = find(strncmp(styles(iref+1:end), 'Heading', 7), 1);
    nrefs = numel(paras) - iref;
    if ~isempty(later)
        fails{end+1} = sprintf('a section follows the references: "%s"', ...
            paras{iref + later}(1:min(60, end)));
    end
    fprintf(['[verify] %d paragraphs, references at %d, %d entries after it, ' ...
             'no later section\n'], numel(paras), iref, nrefs);
end

% With the listings gone, every paragraph is prose and every check applies to
% all of it.
prose = text;

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

% ---- 6. no code listing in the document -----------------------------------
ncode = sum(strcmp(styles, 'SourceCode'));
fprintf('[verify] %d code paragraphs (expected 0)\n', ncode);
if ncode > 0
    fails{end+1} = sprintf(['%d code listing paragraph(s) in the document; ' ...
        'the pipeline belongs in the repository'], ncode);
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
