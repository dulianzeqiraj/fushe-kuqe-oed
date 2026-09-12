function D = fk_load_data(cfg)
%FK_LOAD_DATA  Read the field tables written by prepare_data.py.
%
%  Returns a struct with the well tables, the aquifer boundary polygon and
%  the facies index of every well that carries a pumping-test K.  Nothing is
%  imputed: wells without a measurement keep a NaN.

opts = detectImportOptions(fullfile(cfg.dir_data, 'wells_all.csv'), ...
                           'TextType', 'string');
T = readtable(fullfile(cfg.dir_data, 'wells_all.csv'), opts);

D = struct();
D.well_id  = string(T.Well_ID);
D.x        = T.X_GaussKruger;
D.y        = T.Y_GaussKruger;
D.K        = T.K_m_per_day;            % m/day, NaN where not measured
D.NO3      = T.NO3_mg_per_L;           % mg/L,  NaN where not measured
D.DRASTIC  = T.DRASTIC_Index;          % dimensionless, all 180 wells
D.depth_w  = T.Piezometric_Level_m;    % depth to water, m
D.lith     = string(T.Lithology);
D.inside   = strcmp(string(T.Inside_Aquifer), "YES");

% facies index 1 coarse, 2 medium, 3 fine, 0 if the log gives no lithology
D.facies = zeros(numel(D.well_id), 1);
for f = 1:numel(cfg.facies_names)
    D.facies(D.lith == string(cfg.facies_names{f})) = f;
end

D.has_K   = ~isnan(D.K);
D.has_NO3 = ~isnan(D.NO3);

B = readtable(fullfile(cfg.dir_data, 'boundary.csv'));
D.bx = B.X_GaussKruger;
D.by = B.Y_GaussKruger;

% ---- inventory, printed so that every count in the paper is traceable ----
D.n_total       = numel(D.well_id);
D.n_inside      = sum(D.inside);
D.n_K           = sum(D.has_K);
D.n_NO3         = sum(D.has_NO3);
D.n_inside_K    = sum(D.inside & D.has_K);
D.n_inside_NO3  = sum(D.inside & D.has_NO3);
D.n_boundary    = numel(D.bx);

fprintf('[load] wells total %d, inside %d, with K %d, with NO3 %d\n', ...
        D.n_total, D.n_inside, D.n_K, D.n_NO3);
fprintf('[load] inside with K %d, inside with NO3 %d, boundary vertices %d\n', ...
        D.n_inside_K, D.n_inside_NO3, D.n_boundary);
for f = 1:3
    fprintf('[load] facies %-7s : %d wells with K\n', cfg.facies_short{f}, ...
            sum(D.has_K & D.facies == f));
end
fprintf('[load] K  range %.0f - %.0f m/day, mean %.1f, sd %.1f\n', ...
        min(D.K(D.has_K)), max(D.K(D.has_K)), ...
        mean(D.K(D.has_K)), std(D.K(D.has_K)));
fprintf('[load] NO3 range %.2f - %.2f mg/L, mean %.2f, sd %.2f\n', ...
        min(D.NO3(D.has_NO3)), max(D.NO3(D.has_NO3)), ...
        mean(D.NO3(D.has_NO3)), std(D.NO3(D.has_NO3)));
fprintf('[load] DRASTIC range %g - %g, mean %.2f, sd %.2f (n = %d)\n', ...
        min(D.DRASTIC), max(D.DRASTIC), mean(D.DRASTIC), std(D.DRASTIC), ...
        sum(~isnan(D.DRASTIC)));
end
