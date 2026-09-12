function fk_md2docx(mdtext, template, outfile)
%FK_MD2DOCX  Write a .docx from the manuscript's markdown, in MATLAB alone.
%
%  fk_md2docx(mdtext, template, outfile)
%
%  Handles the subset of markdown the manuscript uses: YAML title and author,
%  ATX headings, paragraphs, bullet and numbered lists, pipe tables, fenced
%  code blocks, bold and code spans, and inline and display mathematics, the
%  last through fk_latex2omml.
%
%  The template supplies every part of the package except word/document.xml:
%  styles, numbering, theme, fonts, section setup.  Only the body is written
%  here, which is the same discipline used elsewhere in this project for
%  editing a Word file without rebuilding it.
%
%  Anything the parser does not recognise raises an error.  Silently dropping
%  a construct would produce a document that looks finished and is not.

W = ['xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main"'];
M = ['xmlns:m="http://schemas.openxmlformats.org/officeDocument/2006/math"'];
R = ['xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"'];

lines = regexp(mdtext, '\r?\n', 'split');
n = numel(lines);
body = {};
i = 1;
after_heading = false;

% ------------------------------------------------------- YAML front matter
if n >= 1 && strcmp(strtrim(lines{1}), '---')
    j = 2;
    ttl = ''; aut = '';
    while j <= n && ~strcmp(strtrim(lines{j}), '---')
        tok = regexp(lines{j}, '^(\w+):\s*"?(.*?)"?\s*$', 'tokens', 'once');
        if ~isempty(tok)
            switch tok{1}
                case 'title',  ttl = tok{2};
                case 'author', aut = tok{2};
            end
        end
        j = j + 1;
    end
    if ~isempty(ttl)
        body{end+1} = local_para('Title', local_inline(ttl));
    end
    if ~isempty(aut)
        body{end+1} = local_para('Author', local_inline(aut));
    end
    i = j + 1;
end

% ------------------------------------------------------------ block loop --
while i <= n
    ln = lines{i};
    s = strtrim(ln);

    if isempty(s)
        i = i + 1;
        continue
    end

    % ---- fenced code ------------------------------------------------------
    % One paragraph for the whole block, with a line break between lines.  A
    % paragraph per line would carry the style's paragraph spacing into every
    % line and more than double the page count.
    if startsWith(s, '```')
        i = i + 1;
        code = {};
        while i <= n && ~startsWith(strtrim(lines{i}), '```')
            code{end+1} = lines{i}; %#ok<AGROW>
            i = i + 1;
        end
        i = i + 1;
        body{end+1} = local_code_block(code); %#ok<AGROW>
        after_heading = false;
        continue
    end

    % ---- pipe table -------------------------------------------------------
    if startsWith(s, '|')
        rows = {};
        while i <= n && startsWith(strtrim(lines{i}), '|')
            rows{end+1} = strtrim(lines{i}); %#ok<AGROW>
            i = i + 1;
        end
        body{end+1} = local_table(rows); %#ok<AGROW>
        after_heading = false;
        continue
    end

    % ---- heading ----------------------------------------------------------
    h = regexp(s, '^(#{1,6})\s+(.*)$', 'tokens', 'once');
    if ~isempty(h)
        lvl = numel(h{1});
        body{end+1} = local_para(sprintf('Heading%d', lvl), ...
                                 local_inline(h{2})); %#ok<AGROW>
        i = i + 1;
        after_heading = true;
        continue
    end

    % ---- display equation on its own line ---------------------------------
    if startsWith(s, '$$') && endsWith(s, '$$') && numel(s) > 4
        tex = s(3:end-2);
        body{end+1} = ['<w:p>' fk_latex2omml(tex, true) '</w:p>']; %#ok<AGROW>
        i = i + 1;
        after_heading = false;
        continue
    end

    % ---- bullet item ------------------------------------------------------
    if startsWith(s, '- ')
        body{end+1} = local_listpara(1001, local_inline(s(3:end))); %#ok<AGROW>
        i = i + 1;
        after_heading = false;
        continue
    end

    % ---- numbered item ----------------------------------------------------
    nm = regexp(s, '^\d+\.\s+(.*)$', 'tokens', 'once');
    if ~isempty(nm)
        body{end+1} = local_listpara(1002, local_inline(nm{1})); %#ok<AGROW>
        i = i + 1;
        after_heading = false;
        continue
    end

    % ---- ordinary paragraph, possibly wrapped over several lines ----------
    buf = {s};
    i = i + 1;
    while i <= n
        s2 = strtrim(lines{i});
        if isempty(s2) || startsWith(s2, '|') || startsWith(s2, '#') || ...
           startsWith(s2, '- ') || startsWith(s2, '```') || ...
           startsWith(s2, '$$') || ~isempty(regexp(s2, '^\d+\.\s', 'once'))
            break
        end
        buf{end+1} = s2; %#ok<AGROW>
        i = i + 1;
    end
    txt = strjoin(buf, ' ');
    if after_heading, sty = 'FirstParagraph'; else, sty = 'BodyText'; end
    body{end+1} = local_para(sty, local_inline(txt)); %#ok<AGROW>
    after_heading = false;
end

% ---------------------------------------------------------- assemble part --
sectPr = local_sectpr(template);
xml = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>' ...
       '<w:document ' W ' ' M ' ' R '><w:body>' ...
       strjoin(body, '') sectPr '</w:body></w:document>'];

fk_zip_replace(template, outfile, 'word/document.xml', xml);
fprintf('[docx] %s written, %d block elements\n', outfile, numel(body));
end

% =========================================================================
function p = local_para(style, runs)
p = ['<w:p><w:pPr><w:pStyle w:val="' style '"/></w:pPr>' runs '</w:p>'];
end

function p = local_listpara(numid, runs)
p = ['<w:p><w:pPr><w:pStyle w:val="Compact"/><w:numPr>' ...
     '<w:ilvl w:val="0"/><w:numId w:val="' num2str(numid) '"/>' ...
     '</w:numPr></w:pPr>' runs '</w:p>'];
end

function p = local_code_block(codelines)
runs = '';
for k = 1:numel(codelines)
    if k > 1
        runs = [runs '<w:r><w:br/></w:r>']; %#ok<AGROW>
    end
    runs = [runs '<w:r><w:rPr><w:rStyle w:val="VerbatimChar"/></w:rPr>' ...
            '<w:t xml:space="preserve">' fk_xmlesc(codelines{k}) ...
            '</w:t></w:r>']; %#ok<AGROW>
end
p = ['<w:p><w:pPr><w:pStyle w:val="SourceCode"/></w:pPr>' runs '</w:p>'];
end

% =========================================================================
function out = local_inline(s)
%LOCAL_INLINE  Split a paragraph into runs: maths, bold, code, plain.
out = '';
rest = s;
while ~isempty(rest)
    % earliest of the four markers
    iM = regexp(rest, '\$[^\$]+\$', 'once');
    iB = regexp(rest, '\*\*[^\*]+\*\*', 'once');
    iC = regexp(rest, '`[^`]+`', 'once');
    iI = regexp(rest, '(?<!\*)\*[^\*]+\*(?!\*)', 'once');
    kinds = {'m', 'b', 'c', 'i'};
    keep = ~cellfun(@isempty, {iM, iB, iC, iI});
    if ~any(keep)
        out = [out local_run(rest, '')];
        return
    end
    cc = [iM iB iC iI];         % empties vanish, so cc lines up with kk
    kk = kinds(keep);
    [pos, w] = min(cc);
    kind = kk{w};

    if pos > 1
        out = [out local_run(rest(1:pos-1), '')];
    end
    rest = rest(pos:end);
    switch kind
        case 'm'
            tk = regexp(rest, '^\$([^\$]+)\$', 'tokens', 'once');
            out = [out fk_latex2omml(tk{1}, false)];
            rest = rest(numel(tk{1}) + 3:end);
        case 'b'
            tk = regexp(rest, '^\*\*([^\*]+)\*\*', 'tokens', 'once');
            out = [out local_run(tk{1}, '<w:b/>')];
            rest = rest(numel(tk{1}) + 5:end);
        case 'c'
            tk = regexp(rest, '^`([^`]+)`', 'tokens', 'once');
            out = [out local_run(tk{1}, '<w:rStyle w:val="VerbatimChar"/>')];
            rest = rest(numel(tk{1}) + 3:end);
        case 'i'
            tk = regexp(rest, '^\*([^\*]+)\*', 'tokens', 'once');
            out = [out local_run(tk{1}, '<w:i/>')];
            rest = rest(numel(tk{1}) + 3:end);
    end
end
end

function r = local_run(txt, rpr)
if isempty(txt), r = ''; return; end
if isempty(rpr)
    r = ['<w:r><w:t xml:space="preserve">' fk_xmlesc(txt) '</w:t></w:r>'];
else
    r = ['<w:r><w:rPr>' rpr '</w:rPr><w:t xml:space="preserve">' ...
         fk_xmlesc(txt) '</w:t></w:r>'];
end
end

% =========================================================================
function t = local_table(rows)
%LOCAL_TABLE  A pipe table.  The second row of a markdown table is the
%             alignment rule and carries no content.
cells = cell(numel(rows), 1);
keep = true(numel(rows), 1);
for k = 1:numel(rows)
    r = rows{k};
    r = regexprep(r, '^\|', '');
    r = regexprep(r, '\|$', '');
    parts = strsplit(r, '|');
    parts = cellfun(@strtrim, parts, 'UniformOutput', false);
    if all(~cellfun(@isempty, regexp(parts, '^:?-{2,}:?$', 'once')))
        keep(k) = false;
    end
    cells{k} = parts;
end
cells = cells(keep);
ncol = max(cellfun(@numel, cells));
w = floor(9360 / ncol);

grid = '';
for c = 1:ncol
    grid = [grid '<w:gridCol w:w="' num2str(w) '"/>']; %#ok<AGROW>
end

trs = '';
for k = 1:numel(cells)
    tcs = '';
    row = cells{k};
    for c = 1:ncol
        if c <= numel(row), txt = row{c}; else, txt = ''; end
        tcs = [tcs '<w:tc><w:tcPr/><w:p><w:pPr><w:pStyle w:val="Compact"/>' ...
               '</w:pPr>' local_inline(txt) '</w:p></w:tc>']; %#ok<AGROW>
    end
    if k == 1
        trs = [trs '<w:tr><w:trPr><w:tblHeader w:val="on"/></w:trPr>' tcs '</w:tr>']; %#ok<AGROW>
    else
        trs = [trs '<w:tr>' tcs '</w:tr>']; %#ok<AGROW>
    end
end

t = ['<w:tbl><w:tblPr><w:tblStyle w:val="Table"/>' ...
     '<w:tblW w:type="pct" w:w="5000"/><w:tblLayout w:type="fixed"/>' ...
     '<w:tblLook w:firstRow="1" w:lastRow="0" w:firstColumn="0" ' ...
     'w:lastColumn="0" w:noHBand="0" w:noVBand="0" w:val="0020"/></w:tblPr>' ...
     '<w:tblGrid>' grid '</w:tblGrid>' trs '</w:tbl>' ...
     '<w:p><w:pPr><w:pStyle w:val="BodyText"/></w:pPr></w:p>'];
end

% =========================================================================
function sp = local_sectpr(template)
%LOCAL_SECTPR  Reuse the template's page setup verbatim.
tmp = tempname;
mkdir(tmp);
c = onCleanup(@() rmdir(tmp, 's'));
unzip(template, tmp);
d = fileread(fullfile(tmp, 'word', 'document.xml'));
m = regexp(d, '<w:sectPr[^>]*>.*?</w:sectPr>', 'match', 'once');
if isempty(m)
    m = regexp(d, '<w:sectPr[^>]*/>', 'match', 'once');
end
if isempty(m)
    error('fk_md2docx:sectPr', 'template has no sectPr');
end
sp = m;
end
