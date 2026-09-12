function fk_prepare_data()
%FK_PREPARE_DATA  Extract the analysis tables from the sources into flat CSVs.
%
%  Inputs, read only and never modified:
%    Fushe_Kuqe_All_Data.xlsx        the Desktop copy, which carries three
%                                    sheets transcribed from Cenameri and
%                                    Beqiraj (2016) besides the well and
%                                    boundary tables
%    Pugh2023_MRVA_44wells_REAL.csv  the USGS Mississippi Alluvial Plain set
%
%  Outputs, into ../data:
%    wells_all.csv      180 wells, every column of the well sheet
%    wells_K.csv        the wells carrying a pumping-test conductivity
%    wells_NO3.csv      the wells carrying a nitrate measurement
%    boundary.csv       the aquifer boundary vertices
%    swi_chemistry.csv  major ions at four wells
%    swi_cl_trend.csv   chloride in 1984, 1999 and 2001
%    swi_facts.csv      flow direction, heads, thicknesses, wedge geometry
%    mrva_44.csv        the USGS set
%    manifest.json      the counts, so a reader can check the extraction
%
%  No value is altered, imputed or rounded: rows are copied or dropped only.

cfg = fk_config();
out = cfg.dir_data;
if ~exist(out, 'dir'), mkdir(out); end

xlsx = 'C:\Users\d_zeq\OneDrive\Desktop\Fushe_Kuqe_All_Data.xlsx';
mrva = ['C:\Users\d_zeq\OneDrive\Desktop\PROJEKTE 2 ARTIKUJ\' ...
        'AKUIFER TE DHENA TE NDRYSHME FUSHE KUQE TE DHENA\' ...
        'Pugh2023_MRVA_44wells_REAL.csv'];
for f = {xlsx, mrva}
    if ~exist(f{1}, 'file')
        error('fk_prepare_data:source', 'missing source %s', f{1});
    end
end

% --------------------------------------------------------------- wells ----
T = readtable(xlsx, 'Sheet', '2_All_Wells', 'VariableNamingRule', 'preserve');
T = T(~ismissing(T.(1)), :);
writetable(T, fullfile(out, 'wells_all.csv'));

hasK = ~ismissing(T.('K_m_per_day'));
hasN = ~ismissing(T.('NO3_mg_per_L'));
writetable(T(hasK, :), fullfile(out, 'wells_K.csv'));
writetable(T(hasN, :), fullfile(out, 'wells_NO3.csv'));

B = readtable(xlsx, 'Sheet', '3_Boundary_Coords', 'VariableNamingRule', 'preserve');
B = B(~ismissing(B.(1)), :);
writetable(B, fullfile(out, 'boundary.csv'));

% ------------------------------------- the three seawater-intrusion sheets --
% Each carries two banner lines and a blank line before the real header, so
% the header row is located rather than assumed.
swi = {'5_Cenameri2016_Chemistry', 'swi_chemistry.csv'; ...
       '6_Cenameri2016_Cl_Trend',  'swi_cl_trend.csv'; ...
       '7_Cenameri2016_SWI_Facts', 'swi_facts.csv'};
sheets = sheetnames(xlsx);
nswi = zeros(size(swi, 1), 1);
for k = 1:size(swi, 1)
    if ~any(strcmp(sheets, swi{k, 1})), continue; end
    C = readcell(xlsx, 'Sheet', swi{k, 1});
    hdr = local_header_row(C);
    S = cell2table(local_clean(C(hdr+1:end, :)), ...
                   'VariableNames', local_names(C(hdr, :)));
    S = S(~local_blank(S), :);
    writetable(S, fullfile(out, swi{k, 2}));
    nswi(k) = height(S);
end

% ---------------------------------------------------------------- MRVA ----
M = readtable(mrva, 'VariableNamingRule', 'preserve');
writetable(M, fullfile(out, 'mrva_44.csv'));

% ------------------------------------------------------------ manifest ----
inside = strcmp(string(T.('Inside_Aquifer')), "YES");
lith = string(T.('Lithology'));
man = struct();
man.source_workbook = xlsx;
man.source_mrva = mrva;
man.n_wells_total = height(T);
man.n_wells_with_K = sum(hasK);
man.n_wells_with_NO3 = sum(hasN);
man.n_wells_inside = sum(inside);
man.n_inside_with_K = sum(inside & hasK);
man.n_inside_with_NO3 = sum(inside & hasN);
man.n_boundary_vertices = height(B);
man.n_mrva = height(M);
man.swi_rows = nswi';
u = unique(lith(hasK));
lc = struct();
for k = 1:numel(u)
    lc.(matlab.lang.makeValidName(char(u(k)))) = sum(hasK & lith == u(k));
end
man.lithology_of_K_wells = lc;

fid = fopen(fullfile(out, 'manifest.json'), 'w');
fwrite(fid, jsonencode(man, 'PrettyPrint', true));
fclose(fid);

fprintf(['[prep] %d wells, %d with K, %d with NO3, %d inside, ' ...
         '%d boundary vertices, %d MRVA\n'], man.n_wells_total, ...
        man.n_wells_with_K, man.n_wells_with_NO3, man.n_wells_inside, ...
        man.n_boundary_vertices, man.n_mrva);
end

% =========================================================================
function h = local_header_row(C)
h = 0;
for k = 2:size(C, 1)
    row = C(k, :);
    filled = sum(cellfun(@(v) ~local_ismissing(v), row));
    if filled >= 2
        h = k;
        return
    end
end
error('fk_prepare_data:header', 'no header row found');
end

function tf = local_ismissing(v)
tf = isa(v, 'missing') || (ischar(v) && isempty(v)) || ...
     (isnumeric(v) && isscalar(v) && isnan(v)) || ...
     (isstring(v) && (ismissing(v) || v == ""));
end

function C = local_clean(C)
for k = 1:numel(C)
    if local_ismissing(C{k}), C{k} = ''; end
    if isstring(C{k}), C{k} = char(C{k}); end
    if isnumeric(C{k}) && isscalar(C{k}), C{k} = num2str(C{k}); end
    if ~ischar(C{k}), C{k} = ''; end
end
end

function nm = local_names(row)
nm = cell(1, numel(row));
for k = 1:numel(row)
    v = row{k};
    if local_ismissing(v), v = sprintf('Var%d', k); end
    if ~ischar(v), v = char(string(v)); end
    nm{k} = matlab.lang.makeValidName(v);
end
nm = matlab.lang.makeUniqueStrings(nm);
end

function tf = local_blank(S)
C = table2cell(S);
tf = all(cellfun(@(v) isempty(v), C), 2);
first = C(:, 1);
tf = tf | cellfun(@(v) ischar(v) && (startsWith(v, 'Note') || ...
                                     startsWith(v, '- ')), first);
end
