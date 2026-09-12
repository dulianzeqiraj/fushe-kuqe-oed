function R = run_all()
%RUN_ALL  Reproduce every number reported in the manuscript.
%
%  Usage:  cd code; matlab -batch "run_all"
%
%  Writes results/results.json, which is the single file the manuscript
%  quotes from, plus results/run_all.mat with the full state.  Any number in
%  the paper that is not in results.json is an error.

t_start = tic;
cfg = fk_config();
if ~exist(cfg.dir_results, 'dir'), mkdir(cfg.dir_results); end

fprintf('\n===== 1. data =====\n');
D = fk_load_data(cfg);

fprintf('\n===== 2. is the nitrate field predictable at all =====\n');
Pr = fk_predictability(cfg, D);

fprintf('\n===== 3. grid and prior =====\n');
G = fk_build_grid(cfg, D);
P = fk_prior_field(cfg, D, G);

fprintf('\n===== 4. forward model and calibration =====\n');
F = fk_forward(cfg, G, P);
Cal = fk_calibrate(cfg, D, G, P, F);

fprintf('\n===== 5. sensitivity and Fisher information =====\n');
S = fk_sensitivity(cfg, G, P, F, Cal);

fprintf('\n===== 6. monitoring network design =====\n');
Des = fk_design(cfg, D, G, P, F, Cal, S);

fprintf('\n===== 7. the same design with a local operator =====\n');
cand = true(G.n_active, 1);
cand(Cal.cell_of) = false;
w = (P.V / mean(P.V)) .^ cfg.gamma_Vdesign;
Loc = oed_local(P.Phi, P.lam, w, cfg.sd_local, cand, cfg.n_new_wells, ...
                G.xa, G.ya, P.vrange);
fprintf(['[local] Fushe-Kuqe, local operator: variance of the weighted mean ' ...
         'down %.1f%% (D) and %.1f%% (goal); median separation %.0f m\n'], ...
        100 * Loc.var_reduction_D, 100 * Loc.var_reduction_goal, ...
        Loc.median_separation_m);

fprintf('\n===== 8. posterior sampling =====\n');
Par = fk_particles(cfg, G, P, F, Cal, S);

fprintf('\n===== 9. independent replication on the USGS MRVA set =====\n');
M = mrva_replicate(cfg);

fprintf('\n===== 10. sensitivity of the conclusions =====\n');
Sw = fk_sweep(cfg, D, G, P, F, Cal, S, Des);

R = struct('cfg', cfg, 'D', D, 'Pr', Pr, 'G', G, 'P', P, 'F', F, ...
           'Cal', Cal, 'S', S, 'Des', Des, 'Loc', Loc, 'Par', Par, ...
           'M', M, 'Sw', Sw);
R.runtime_s = toc(t_start);
R.matlab_version = version;
R.run_date = datestr(now, 'yyyy-mm-dd HH:MM:SS');

save(fullfile(cfg.dir_results, 'run_all.mat'), '-struct', 'R', '-v7.3');
fk_write_results(cfg, R);
fprintf('\nfinished in %.0f s\n', R.runtime_s);
end
