function xml = fk_latex2omml(tex, display)
%FK_LATEX2OMML  Convert the LaTeX subset used in this manuscript to OMML.
%
%  xml = fk_latex2omml(tex)           inline equation
%  xml = fk_latex2omml(tex, true)     display equation, wrapped in m:oMathPara
%
%  This is not a general LaTeX engine.  It covers exactly the constructs the
%  manuscript uses, which fk_check_omml.m enumerates from the source and
%  asserts are all handled; anything outside that set raises an error rather
%  than being silently dropped, because a silently dropped symbol in an
%  equation is the kind of error that survives to print.
%
%  Word stores equations as Office MathML.  Writing it directly keeps the
%  equations editable in Word, which rendering them as images would not.

if nargin < 2, display = false; end

toks = local_tokenise(tex);
[body, pos] = local_parse(toks, 1, numel(toks));
if pos <= numel(toks)
    error('fk_latex2omml:trailing', 'unconsumed input in "%s"', tex);
end

if display
    xml = ['<m:oMathPara><m:oMath>' body '</m:oMath></m:oMathPara>'];
else
    xml = ['<m:oMath>' body '</m:oMath>'];
end
end

% =========================================================================
function t = local_tokenise(s)
%LOCAL_TOKENISE  Split into command, brace, script and character tokens.
t = {};
i = 1;
n = numel(s);
while i <= n
    c = s(i);
    if c == '\'
        if i == n, error('fk_latex2omml:dangling', 'trailing backslash'); end
        j = i + 1;
        if isletter(s(j))
            while j <= n && isletter(s(j)), j = j + 1; end
            t{end+1} = s(i:j-1); %#ok<AGROW>
            i = j;
        else
            t{end+1} = s(i:j); %#ok<AGROW>
            i = j + 1;
        end
    elseif any(c == '{}^_')
        t{end+1} = c; %#ok<AGROW>
        i = i + 1;
    elseif c == ' '
        i = i + 1;                      % spacing is Word's business
    else
        t{end+1} = c; %#ok<AGROW>
        i = i + 1;
    end
end
end

% =========================================================================
function [out, i] = local_parse(t, i, stop)
%LOCAL_PARSE  Parse tokens i..stop into OMML, collapsing adjacent characters
%             into single runs so Word does not get one run per letter.
out = '';
buf = '';
    function flush()
        if ~isempty(buf)
            out = [out local_run(buf)];
            buf = '';
        end
    end

while i <= stop
    tk = t{i};
    if strcmp(tk, '}')
        break
    elseif strcmp(tk, '{')
        [inner, i] = local_parse(t, i + 1, stop);
        if i > stop || ~strcmp(t{i}, '}')
            error('fk_latex2omml:brace', 'unbalanced brace');
        end
        i = i + 1;
        [inner, i] = local_scripts(t, i, stop, inner);
        flush();
        out = [out inner];
    elseif strcmp(tk, '^') || strcmp(tk, '_')
        % a script attached to the characters accumulated so far
        if isempty(buf)
            error('fk_latex2omml:script', 'script with no base');
        end
        base = local_run(buf(end));
        buf(end) = [];
        flush();
        [sc, i] = local_scripts(t, i, stop, base);
        out = [out sc];
    else
        [node, i, isrun] = local_atom(t, i, stop);
        if isrun
            buf = [buf node];
        else
            [node, i] = local_scripts(t, i, stop, node);
            flush();
            out = [out node];
        end
    end
end
flush();
end

% =========================================================================
function [out, i] = local_scripts(t, i, stop, base)
%LOCAL_SCRIPTS  Attach any ^ and _ that follow a completed base.
sup = ''; sub = '';
while i <= stop && (strcmp(t{i}, '^') || strcmp(t{i}, '_'))
    kind = t{i};
    [arg, i] = local_argument(t, i + 1, stop);
    if strcmp(kind, '^'), sup = arg; else, sub = arg; end
end
if ~isempty(sup) && ~isempty(sub)
    out = ['<m:sSubSup><m:e>' base '</m:e><m:sub>' sub ...
           '</m:sub><m:sup>' sup '</m:sup></m:sSubSup>'];
elseif ~isempty(sup)
    out = ['<m:sSup><m:e>' base '</m:e><m:sup>' sup '</m:sup></m:sSup>'];
elseif ~isempty(sub)
    out = ['<m:sSub><m:e>' base '</m:e><m:sub>' sub '</m:sub></m:sSub>'];
else
    out = base;
end
end

% =========================================================================
function [out, i] = local_argument(t, i, stop)
%LOCAL_ARGUMENT  One braced group, or one atom.
if i > stop, error('fk_latex2omml:arg', 'missing argument'); end
if strcmp(t{i}, '{')
    [out, i] = local_parse(t, i + 1, stop);
    if i > stop || ~strcmp(t{i}, '}')
        error('fk_latex2omml:brace', 'unbalanced brace in argument');
    end
    i = i + 1;
else
    [node, i, isrun] = local_atom(t, i, stop);
    if isrun, out = local_run(node); else, out = node; end
end
end

% =========================================================================
function [out, i, isrun] = local_atom(t, i, stop)
%LOCAL_ATOM  One token, expanded.  isrun true means "plain text, may be
%            merged with neighbours into one run".
tk = t{i};
isrun = false;

% ---- structures taking arguments -----------------------------------------
switch tk
    case '\frac'
        [num, i] = local_argument(t, i + 1, stop);
        [den, i] = local_argument(t, i, stop);
        out = ['<m:f><m:fPr><m:type m:val="bar"/></m:fPr>' ...
               '<m:num>' num '</m:num><m:den>' den '</m:den></m:f>'];
        return
    case '\sqrt'
        [e, i] = local_argument(t, i + 1, stop);
        out = ['<m:rad><m:radPr><m:degHide m:val="1"/></m:radPr>' ...
               '<m:deg/><m:e>' e '</m:e></m:rad>'];
        return
    case {'\bar', '\overline'}
        [e, i] = local_argument(t, i + 1, stop);
        out = ['<m:bar><m:barPr><m:pos m:val="top"/></m:barPr>' ...
               '<m:e>' e '</m:e></m:bar>'];
        return
    case {'\text', '\mathcal', '\mathrm'}
        [e, i] = local_argument(t, i + 1, stop);
        out = e;                       % upright is Word's default for m:t
        return
    case '\operatorname'
        % An upright run, not an m:func.  A function with an empty base is
        % valid OMML but Word reserves space for the missing argument and
        % leaves a gap before the bracket that follows.
        j = i + 1;
        if j <= stop && strcmp(t{j}, '{')
            k = j + 1; txt = '';
            while k <= stop && ~strcmp(t{k}, '}')
                txt = [txt t{k}]; %#ok<AGROW>
                k = k + 1;
            end
            i = k + 1;
        else
            txt = t{j}; i = j + 1;
        end
        out = local_upright(txt);
        return
    case '\underbrace'
        [e, i] = local_argument(t, i + 1, stop);
        if i <= stop && strcmp(t{i}, '_')
            [lab, i] = local_argument(t, i + 1, stop);
        else
            lab = '';
        end
        grp = ['<m:groupChr><m:groupChrPr><m:chr m:val="&#9183;"/>' ...
               '<m:pos m:val="bot"/><m:vertJc m:val="top"/></m:groupChrPr>' ...
               '<m:e>' e '</m:e></m:groupChr>'];
        if isempty(lab)
            out = grp;
        else
            out = ['<m:limLow><m:e>' grp '</m:e><m:lim>' lab '</m:lim></m:limLow>'];
        end
        return
    case '\sum'
        sub = ''; sup = '';
        while i + 1 <= stop && (strcmp(t{i+1}, '_') || strcmp(t{i+1}, '^'))
            kind = t{i+1};
            [arg, i2] = local_argument(t, i + 2, stop);
            if strcmp(kind, '_'), sub = arg; else, sup = arg; end
            i = i2 - 1;
        end
        i = i + 1;
        if isempty(sup)
            % A sum with only a lower index is set as the symbol with a
            % subscript.  An n-ary with an empty base would be correct OMML
            % but Word reserves space for the missing summand and leaves a
            % visible gap before the terms that follow.
            base = local_run(char(8721));
            if isempty(sub)
                out = base;
            else
                out = ['<m:sSub><m:e>' base '</m:e><m:sub>' sub ...
                       '</m:sub></m:sSub>'];
            end
        else
            pr = ['<m:naryPr><m:chr m:val="&#8721;"/>' ...
                  '<m:limLoc m:val="undOvr"/>'];
            if isempty(sub), pr = [pr '<m:subHide m:val="1"/>']; end
            pr = [pr '</m:naryPr>'];
            out = ['<m:nary>' pr '<m:sub>' sub '</m:sub><m:sup>' sup ...
                   '</m:sup><m:e></m:e></m:nary>'];
        end
        return
    case '\left'
        i = i + 1;
        if i > stop, error('fk_latex2omml:left', '\\left with no delimiter'); end
        opench = local_delim(t{i});
        i = i + 1;
        depth = 0;
        j = i;
        while j <= stop
            if strcmp(t{j}, '\left'), depth = depth + 1; end
            if strcmp(t{j}, '\right')
                if depth == 0, break; end
                depth = depth - 1;
            end
            j = j + 1;
        end
        if j > stop, error('fk_latex2omml:right', 'missing \\right'); end
        [inner, ~] = local_parse(t, i, j - 1);
        closech = local_delim(t{j + 1});
        i = j + 2;
        out = ['<m:d><m:dPr><m:begChr m:val="' opench '"/><m:endChr m:val="' ...
               closech '"/></m:dPr><m:e>' inner '</m:e></m:d>'];
        return
end

% ---- functions rendered upright ------------------------------------------
fn = {'\log', '\ln', '\exp', '\det', '\max', '\min', '\inf', '\sup'};
if any(strcmp(tk, fn))
    name = tk(2:end);
    i = i + 1;
    sub = '';
    if i <= stop && strcmp(t{i}, '_')
        [sub, i] = local_argument(t, i + 1, stop);
    end
    out = local_upright([name char(8201)]);   % thin space after the name
    if ~isempty(sub)
        out = ['<m:sSub><m:e>' local_upright(name) '</m:e><m:sub>' sub ...
               '</m:sub></m:sSub>' local_upright(char(8201))];
    end
    return
end

% ---- single symbols -------------------------------------------------------
% found and sym are separate because several commands map to the empty
% string: \! is negative thin space, which Word handles by itself.
[found, sym] = local_symbol(tk);
if found
    i = i + 1;
    out = sym;
    isrun = true;
    return
end

if numel(tk) == 1
    i = i + 1;
    out = tk;
    isrun = true;
    return
end

error('fk_latex2omml:unknown', 'unhandled LaTeX token "%s"', tk);
end

% =========================================================================
function d = local_delim(tk)
switch tk
    case {'(', ')'}, d = tk;
    case {'[', ']'}, d = tk;
    case {'\{'},     d = '{';
    case {'\}'},     d = '}';
    case {'|'},      d = '|';
    case {'.'},      d = '';
    otherwise, error('fk_latex2omml:delim', 'unsupported delimiter "%s"', tk);
end
end

% =========================================================================
function [found, s] = local_symbol(tk)
%LOCAL_SYMBOL  Unicode for the commands the manuscript uses.  found is false
%              for anything else, so local_atom can raise a clear error; it is
%              separate from s because some commands map to the empty string.
persistent map
if isempty(map)
    map = containers.Map('KeyType', 'char', 'ValueType', 'char');
    g = { '\alpha',char(945); '\beta',char(946); '\gamma',char(947); ...
          '\Gamma',char(915); '\delta',char(948); '\Delta',char(916); ...
          '\epsilon',char(949); '\eta',char(951); '\theta',char(952); ...
          '\Theta',char(920); '\kappa',char(954); '\lambda',char(955); ...
          '\Lambda',char(923); '\mu',char(956); '\nu',char(957); ...
          '\xi',char(958); '\Xi',char(926); '\pi',char(960); ...
          '\Pi',char(928); '\rho',char(961); '\sigma',char(963); ...
          '\Sigma',char(931); '\tau',char(964); '\phi',char(966); ...
          '\Phi',char(934); '\chi',char(967); '\psi',char(968); ...
          '\Psi',char(936); '\omega',char(969); '\Omega',char(937); ...
          '\partial',char(8706); '\ell',char(8467); '\top',char(8868); ...
          '\in',char(8712); '\pm',char(177); '\sim',char(8764); ...
          '\mapsto',char(8614); '\square',char(9633); '\times',char(215); ...
          '\cdot',char(183); '\leq',char(8804); '\geq',char(8805); ...
          '\neq',char(8800); '\approx',char(8776); '\infty',char(8734); ...
          '\{','{'; '\}','}'; '\,',' '; '\;',' '; '\ ',' '; ...
          '\!',''; '\qquad','    '; '\quad','  '; '\%','%' };
    for k = 1:size(g, 1)
        map(g{k, 1}) = g{k, 2};
    end
end
found = isKey(map, tk);
if found
    s = map(tk);
else
    s = '';
end
end

% =========================================================================
function r = local_run(txt)
if isempty(txt)
    r = '';
else
    r = ['<m:r><m:t xml:space="preserve">' fk_xmlesc(txt) '</m:t></m:r>'];
end
end

% =========================================================================
function r = local_upright(txt)
%LOCAL_UPRIGHT  A run set in roman, for operator and function names, which
%               Word would otherwise italicise as if they were variables.
r = ['<m:r><m:rPr><m:nor/></m:rPr><m:t xml:space="preserve">' ...
     fk_xmlesc(txt) '</m:t></m:r>'];
end
