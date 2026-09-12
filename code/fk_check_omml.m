function fk_check_omml(mdfile)
%FK_CHECK_OMML  Convert every equation in the manuscript and report failures.
%
%  The converter raises on anything it does not recognise rather than dropping
%  it, so running it over the whole manuscript is the test that it covers the
%  source.  Called by fk_build_docx before it writes anything.

if nargin < 1
    cfg = fk_config();
    mdfile = fullfile(cfg.dir_code, '..', 'manuscript', 'manuscript_filled.md');
end
s = fileread(mdfile);

disp_eq = regexp(s, '\$\$(.*?)\$\$', 'tokens');
rest = regexprep(s, '\$\$.*?\$\$', '');
inl_eq = regexp(rest, '(?<!\$)\$([^\$\n]+?)\$(?!\$)', 'tokens');

nd = numel(disp_eq);  ni = numel(inl_eq);
fprintf('[omml] %d display and %d inline equations\n', nd, ni);

bad = 0;
for k = 1:nd
    try
        x = fk_latex2omml(disp_eq{k}{1}, true);
        assert(contains(x, 'm:oMathPara'));
    catch ME
        bad = bad + 1;
        fprintf('[omml] DISPLAY %d failed: %s\n   %s\n', k, ME.message, ...
                strtrim(disp_eq{k}{1}));
    end
end
for k = 1:ni
    try
        fk_latex2omml(inl_eq{k}{1}, false);
    catch ME
        bad = bad + 1;
        fprintf('[omml] INLINE %d failed: %s\n   %s\n', k, ME.message, ...
                strtrim(inl_eq{k}{1}));
    end
end

if bad > 0
    error('fk_check_omml:failed', '%d of %d equations did not convert', ...
          bad, nd + ni);
end
fprintf('[omml] all %d equations convert\n', nd + ni);
end
