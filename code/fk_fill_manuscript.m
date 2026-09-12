function fk_fill_manuscript(github_url, zenodo_doi)
%FK_FILL_MANUSCRIPT  Put the run's numbers into the manuscript.
%
%  fk_fill_manuscript()
%  fk_fill_manuscript(github_url, zenodo_doi)
%
%  Substitutes every {{placeholder}} in manuscript/manuscript.md with a value
%  taken from results/results.json and writes manuscript/manuscript_filled.md.
%  It refuses to write anything if a placeholder has no value behind it, and
%  reports values that were computed and never cited, which catches a claim
%  dropped from the text but not from the code.
%
%  This is what keeps the manuscript honest: a number reaches the text only by
%  coming out of a run.
%
%  The output is also checked for control characters. That check is here
%  because one crept into the source on 12 September 2026 and turned
%  $\alpha_L$ into $\x07lpha_L$; the converter dropped the byte silently and
%  the equation reached a delivered document as the word "lpha".

cfg = fk_config();
man = fullfile(cfg.dir_code, '..', 'manuscript');
src = fullfile(man, 'manuscript.md');
dst = fullfile(man, 'manuscript_filled.md');
res = fullfile(cfg.dir_results, 'results.json');

if nargin < 1 || isempty(github_url)
    github_url = 'https://github.com/dulianzeqiraj/fushe-kuqe-oed';
end
if nargin < 2 || isempty(zenodo_doi)
    zenodo_doi = 'https://doi.org/10.5281/zenodo.PLACEHOLDER';
end

r = jsondecode(fileread(res));
m = local_values(r, github_url, zenodo_doi);

txt = fileread(src);
keys = fieldnames(m);
used = false(numel(keys), 1);
for k = 1:numel(keys)
    tag = ['{{' keys{k} '}}'];
    if contains(txt, tag)
        used(k) = true;
        txt = strrep(txt, tag, m.(keys{k}));
    end
end

left = regexp(txt, '\{\{(\w+)\}\}', 'tokens');
if ~isempty(left)
    names = unique(cellfun(@(c) c{1}, left, 'UniformOutput', false));
    error('fk_fill_manuscript:unresolved', ...
          'no value for: %s', strjoin(names, ', '));
end

ctrl = double(txt) < 32 & ~ismember(double(txt), [9 10 13]);
if any(ctrl)
    j = find(ctrl, 1);
    error('fk_fill_manuscript:control', ...
          'control character U+%04X at offset %d, context "%s"', ...
          double(txt(j)), j, txt(max(1, j-40):min(numel(txt), j+20)));
end

fid = fopen(dst, 'w', 'n', 'UTF-8');
fwrite(fid, unicode2native(txt, 'UTF-8'));
fclose(fid);

fprintf('[fill] %s written, %d placeholders resolved\n', dst, sum(used));
if any(~used)
    fprintf('[fill] computed but unused: %s\n', ...
            strjoin(keys(~used)', ', '));
end
end

% =========================================================================
function m = local_values(r, github_url, zenodo_doi)
d = r.data;  p = r.predictability;  g = r.grid;  pr = r.prior;
fw = r.forward;  cal = r.calibration;  fi = r.fisher;  des = r.design;
loc = r.local;  sam = r.sampling;  mr = r.mrva;  sw = r.sweep;
m = struct();
f = @(fmt, v) sprintf(fmt, v);
pct = @(x) [sprintf('%.1f', 100*x) ' per cent'];
rng2 = @(v, fmt) [sprintf(fmt, min(v)) ' to ' sprintf(fmt, max(v))];

% ---------------------------------------------------------------- data ---
m.n_total = f('%d', d.n_wells_total);
m.n_inside = f('%d', d.n_wells_inside);
m.n_K = f('%d', d.n_with_K);
m.n_NO3 = f('%d', d.n_with_NO3);
m.n_inside_K = f('%d', d.n_inside_with_K);
m.n_inside_NO3 = f('%d', d.n_inside_with_NO3);
m.n_boundary = f('%d', d.n_boundary_vertices);
m.K_min = f('%.0f', d.K_min);   m.K_max = f('%.0f', d.K_max);
m.K_mean = f('%.1f', d.K_mean); m.K_sd = f('%.1f', d.K_sd);
m.NO3_min = f('%.2f', d.NO3_min);   m.NO3_max = f('%.2f', d.NO3_max);
m.NO3_mean = f('%.2f', d.NO3_mean); m.NO3_sd = f('%.2f', d.NO3_sd);
m.DR_min = f('%.0f', d.DRASTIC_min); m.DR_max = f('%.0f', d.DRASTIC_max);
m.DR_mean = f('%.1f', d.DRASTIC_mean);
m.facies_c = f('%d', d.facies_counts(1));
m.facies_m = f('%d', d.facies_counts(2));
m.facies_f = f('%d', d.facies_counts(3));
m.eta2 = f('%.3f', d.facies_K_eta2);
m.Fstat = f('%.2f', d.facies_K_F);
m.pval = f('%.2f', d.facies_K_p);
m.df1 = '2';
m.df2 = f('%d', d.n_with_K - 3);

% ------------------------------------------------------ predictability ---
m.best_pred_R2 = f('%+.2f', p.best_R2);
m.krig_min = f('%+.2f', min(p.kriging_R2));
m.krig_max = f('%+.2f', max(p.kriging_R2));
m.idw_min = f('%+.2f', min(p.idw_R2));
m.idw_max = f('%+.2f', max(p.idw_R2));
lin = local_lookup(p.linear_names, p.linear_R2);
m.lin_xy = f('%+.2f', lin('X,Y'));
m.lin_drastic = f('%+.2f', lin('DRASTIC'));
m.corr_north = f('%+.2f', p.corr_northing);
m.corr_along = f('%+.2f', p.corr_alongflow);
m.corr_drastic = f('%+.2f', p.corr_drastic);
m.corr_K = f('%+.2f', p.corr_K);

% ------------------------------------------------------ grid and prior ---
m.nx = f('%d', g.nx);  m.ny = f('%d', g.ny);
m.dx = f('%g', g.dx_m);
m.n_active = f('%d', g.n_active);
m.area_grid = f('%.1f', g.area_km2);
m.area_poly = f('%.1f', g.area_polygon_km2);
m.ne_mean = f('%.3f', pr.network_mean);
m.ne_sd = f('%.3f', pr.network_sd);
m.sd_iv = f('%.4f', pr.sd_independent);
m.sd_c = f('%.4f', pr.sd_facies_correlated(1));
m.sd_m = f('%.4f', pr.sd_facies_correlated(2));
m.sd_f = f('%.4f', pr.sd_facies_correlated(3));
m.sd_used = f('%.4f', pr.sd_used);
m.sd_inflation = f('%.1f', pr.sd_inflation_pct);
m.nugget = f('%.4f', pr.variogram_nugget);
m.sill = f('%.4f', pr.variogram_sill);
m.vrange = f('%.0f', pr.variogram_range_m);
m.n_kl = f('%d', pr.n_kl);
m.kl_var = f('%.0f', 100*pr.kl_variance_captured);

% ------------------------------------------------------------- forward ---
m.gradient = f('%.2e', fw.gradient);
m.q0 = f('%.4f', fw.darcy_flux_m_per_day);
m.gam_min = f('%.2f', fw.Gamma_min);
m.gam_max = f('%.2f', fw.Gamma_max);

% --------------------------------------------------------- calibration ---
m.T0_yr = f('%.1f', cal.T0_years);
m.cal_rmse = f('%.2f', cal.rmse_mg_per_L);
m.cal_R2 = f('%.2f', cal.R2);
m.cal_loo_rmse = f('%.2f', cal.loo_rmse_mg_per_L);
m.cal_loo_R2 = f('%+.2f', cal.loo_R2);
m.drastic_contrib = f('%.3f', cal.R2 - cal.R2_without_DRASTIC);
m.null_R2 = f('%.2f', cal.null_best_R2);
m.null_kappa = f('%.2f', cal.null_best_kappa);
m.frac_neg = [sprintf('%.0f', 100*cal.frac_source_negative_before_clipping) ' per cent'];

% -------------------------------------------------- Fisher information ---
m.tr_obs = f('%.3e', fi.trace_observations);
m.tr_prior = f('%.3e', fi.trace_prior);
m.frac_obs_sci = f('%.1e', fi.observation_share);
m.dofs = f('%.4f', fi.dofs);
m.n_modes = f('%d', fi.n_modes);
m.info_gain = f('%.3f', fi.information_gain_nats);
m.sd_red_pct = f('%.2f', 100*fi.sd_reduction_fraction);
m.sd_prior_mean = f('%.5f', fi.sd_prior_mean);
m.sd_post_mean = f('%.5f', fi.sd_post_mean);
m.fd_check = f('%.1e', fi.finite_difference_check);
m.fd_order = f('%.2f', fi.fd_convergence_order);

% -------------------------------------------------------------- design ---
m.horizon_frac = f('%.0f', 100*des.horizon_fraction);
m.n_new = f('%d', des.n_new_wells);
m.sd_local = f('%.3f', des.sd_local);
m.tau_mean_yr = f('%.1f', des.tau_weighted_mean_days/365.25);
m.sd_goal_prior = f('%.1f', des.sd_goal_prior_days);
m.sd_goal_D = f('%.1f', des.sd_goal_after_D_days);
m.sd_goal_G = f('%.1f', des.sd_goal_after_goal_days);
m.sd_goal_D_hor = f('%.1f', des.sd_goal_after_D_horizon_days);
m.sd_goal_G_hor = f('%.1f', des.sd_goal_after_goal_horizon_days);
m.sd_goal_red_pct = f('%.1f', 100*(1 - des.sd_goal_after_goal_days/des.sd_goal_prior_days));
m.median_sep_m = f('%.0f', des.median_separation_m);
m.V_goal = f('%.1f', mean(des.goal_DRASTIC));
m.V_D = f('%.1f', mean(des.D_DRASTIC));
m.V_domain = f('%.1f', des.domain_mean_DRASTIC);
adv = (des.sd_goal_after_D_days^2 - des.sd_goal_after_goal_days^2) / ...
      des.sd_goal_after_D_days^2;
m.goal_adv_pct = pct(adv);
m.horizon_med = f('%.0f', median(des.information_horizon_m));
m.horizon_min = f('%.0f', min(des.information_horizon_m));
m.horizon_max = f('%.0f', max(des.information_horizon_m));
m.horizon_over_range = f('%.1f', des.information_horizon_median_m / pr.variogram_range_m);
m.corr_constrained = f('%.3f', des.information_correlation_max);
m.corr_horizon = f('%.3f', des.information_correlation_horizon_max);
m.corr_free = f('%.3f', des.information_correlation_unconstrained_max);
m.free_span = f('%.0f', des.unconstrained_span_m);
m.min_sep = f('%.0f', des.min_separation_m);
m.R_equiv = f('%.2f', des.R_equivalence);
m.R_inverse = f('%.0f', 1/des.R_equivalence);
m.sd_pump_ne = f('%.4f', des.sd_pump_test_on_ne);
m.gain_conc = f('%.3g', des.gain_conc_best);
m.gain_pump = f('%.3g', des.gain_pump_best);

% ------------------------------------------------------ local operator ---
m.loc_D_pct = pct(loc.var_reduction_D);
m.loc_goal_pct = pct(loc.var_reduction_goal);
m.loc_ratio = f('%.1f', loc.var_reduction_goal / loc.var_reduction_D);
m.loc_sep_m = f('%.0f', loc.median_separation_m);
m.loc_sep_km = f('%.1f', loc.median_separation_m/1000);

% ------------------------------------------------------------ sampling ---
m.gn_iters = f('%d', sam.gn_iterations);
m.cond_H = f('%.1f', sam.condition_number);
m.iters_nat = f('%d', sam.iters_natural);
m.iters_euc = f('%d', sam.iters_euclidean);
m.sqrt_kappa = f('%.1f', sam.sqrt_condition);
m.M_imp = f('%d', sam.M_implicit);
m.ess_imp = f('%.3f', sam.ess_implicit);
m.M_sir = f('%d', sam.M_sir);
m.ess_sir = f('%.3f', sam.ess_sir);

% ---------------------------------------------------------------- MRVA ---
m.mrva_n = f('%d', mr.n_wells);
m.mrva_K_min = f('%.2f', mr.K_min);
m.mrva_K_max = f('%.1f', mr.K_max);
m.mrva_K_mean = f('%.1f', mr.K_mean);
m.mrva_V_min = f('%.0f', mr.V_min);
m.mrva_V_max = f('%.0f', mr.V_max);
m.mrva_mismatch = f('%d', mr.drastic_mismatch);
m.mrva_sep_km = f('%.1f', mr.median_separation_m/1000);
m.mrva_D_pct = pct(mr.var_reduction_D);
m.mrva_goal_pct = pct(mr.var_reduction_goal);
m.mrva_V_goal = f('%.1f', mr.mean_DRASTIC_goal);
m.mrva_V_D = f('%.1f', mr.mean_DRASTIC_D);

% --------------------------------------------------------------- sweep ---
m.sw_frac_range = rng2(sw.frac_obs_vs_sigma, '%.1e');
m.sw_R_sigma_range = rng2(sw.R_vs_sigma, '%.2f');
m.sw_R_sigmaK_range = rng2(sw.R_vs_sigmaK, '%.2f');
m.sw_sep_range = rng2(sw.sep_vs_gamma, '%.0f');
m.sw_R_alpha_range = rng2(sw.R_vs_alphaL, '%.2f');

% -------------------------------------------------------------- tables ---
m.TABLE_INVENTORY = local_tbl_inventory(d, g, pr);
m.TABLE_LOO = local_tbl_loo(p, cal);
m.TABLE_DESIGN = local_tbl_design(des, loc, mr);
m.TABLE_SWEEP = local_tbl_sweep(sw);

m.GITHUB_URL = github_url;
m.ZENODO_DOI = zenodo_doi;
end

% =========================================================================
function h = local_lookup(names, vals)
h = containers.Map('KeyType', 'char', 'ValueType', 'double');
for k = 1:numel(names)
    nm = names{k};
    if isstring(nm) || iscell(nm), nm = char(nm); end
    h(nm) = vals(k);
end
end

% =========================================================================
function t = local_tbl_inventory(d, g, pr)
rows = {
 'Monitoring points with a DRASTIC index', sprintf('%d', d.n_wells_total)
 'Of which inside the digitised boundary', sprintf('%d', d.n_wells_inside)
 'With a pumping-test hydraulic conductivity', sprintf('%d', d.n_with_K)
 'With a nitrate concentration', sprintf('%d', d.n_with_NO3)
 'Inside the boundary, with conductivity', sprintf('%d', d.n_inside_with_K)
 'Inside the boundary, with nitrate', sprintf('%d', d.n_inside_with_NO3)
 'Aquifer boundary vertices', sprintf('%d', d.n_boundary_vertices)
 ['Digitised area (km' char(178) ')'], sprintf('%.1f', g.area_polygon_km2)
 ['Hydraulic conductivity (m d' char(8315) char(185) ')'], ...
     sprintf('%.0f to %.0f, mean %.1f, sd %.1f', d.K_min, d.K_max, d.K_mean, d.K_sd)
 ['Nitrate (mg L' char(8315) char(185) ')'], ...
     sprintf('%.2f to %.2f, mean %.2f, sd %.2f', d.NO3_min, d.NO3_max, d.NO3_mean, d.NO3_sd)
 'DRASTIC index', ...
     sprintf('%.0f to %.0f, mean %.1f, sd %.1f', d.DRASTIC_min, d.DRASTIC_max, d.DRASTIC_mean, d.DRASTIC_sd)
 'Lithological classes of the K wells (coarse/medium/fine)', ...
     sprintf('%d / %d / %d', d.facies_counts(1), d.facies_counts(2), d.facies_counts(3))
 'ANOVA of K on the lithological label', ...
     sprintf('%s%s = %.3f, F(2,%d) = %.2f, p = %.2f', char(951), char(178), ...
             d.facies_K_eta2, d.n_with_K-3, d.facies_K_F, d.facies_K_p)
 ['Exponential variogram of log' char(8321) char(8320) 'K'], ...
     sprintf('nugget %.4f, partial sill %.4f, range %.0f m', ...
             pr.variogram_nugget, pr.variogram_sill, pr.variogram_range_m)
 };
t = local_mdtable({'Quantity', 'Value'}, rows);
end

function t = local_tbl_loo(p, cal)
rows = {};
for k = 1:numel(p.idw_power)
    rows(end+1, :) = {sprintf('Inverse distance, power %d', p.idw_power(k)), ...
                      sprintf('%+.3f', p.idw_R2(k))}; %#ok<AGROW>
end
for k = 1:numel(p.kriging_range_m)
    rows(end+1, :) = {sprintf('Ordinary kriging, range %g km', p.kriging_range_m(k)/1000), ...
                      sprintf('%+.3f', p.kriging_R2(k))}; %#ok<AGROW>
end
for k = 1:numel(p.linear_names)
    nm = p.linear_names{k};
    if isstring(nm), nm = char(nm); end
    rows(end+1, :) = {sprintf('Linear regression on %s', nm), ...
                      sprintf('%+.3f', p.linear_R2(k))}; %#ok<AGROW>
end
rows(end+1, :) = {'Transient transport model', sprintf('%+.3f', cal.loo_R2)};
rows(end+1, :) = {'Network mean (baseline)', '0.000'};
t = local_mdtable({'Method', 'Leave-one-out R2'}, rows);
end

function t = local_tbl_design(des, loc, mr)
rows = {
 'Transport observations, Fushe-Kuqe', ...
   sprintf('%.1f', 100*(1 - (des.sd_goal_after_D_days/des.sd_goal_prior_days)^2)), ...
   sprintf('%.1f', 100*(1 - (des.sd_goal_after_goal_days/des.sd_goal_prior_days)^2)), ...
   sprintf('%.0f', des.median_separation_m)
 'Characterisation points, Fushe-Kuqe', ...
   sprintf('%.1f', 100*loc.var_reduction_D), ...
   sprintf('%.1f', 100*loc.var_reduction_goal), ...
   sprintf('%.0f', loc.median_separation_m)
 'Characterisation points, MRVA', ...
   sprintf('%.1f', 100*mr.var_reduction_D), ...
   sprintf('%.1f', 100*mr.var_reduction_goal), ...
   sprintf('%.0f', mr.median_separation_m)
 };
t = local_mdtable({'Observation type and site', 'D-optimal (%)', ...
                   'Goal-oriented (%)', 'Median separation (m)'}, rows);
end

function t = local_tbl_sweep(sw)
rows = {
 ['Nitrate measurement sd (mg L' char(8315) char(185) ')'], ...
   sprintf('%g to %g', min(sw.sigma_obs), max(sw.sigma_obs)), ...
   sprintf('%.3f to %.3f', min(sw.R_vs_sigma), max(sw.R_vs_sigma)), ...
   sprintf('observation share %.1e to %.1e', min(sw.frac_obs_vs_sigma), max(sw.frac_obs_vs_sigma))
 'Pumping-test relative error', ...
   sprintf('%g to %g', min(sw.sigma_K_rel), max(sw.sigma_K_rel)), ...
   sprintf('%.3f to %.3f', min(sw.R_vs_sigmaK), max(sw.R_vs_sigmaK)), '-'
 'Longitudinal dispersivity (m)', ...
   sprintf('%g to %g', min(sw.alphaL), max(sw.alphaL)), ...
   sprintf('%.3f to %.3f', min(sw.R_vs_alphaL), max(sw.R_vs_alphaL)), ...
   sprintf('median separation %.0f to %.0f m', min(sw.sep_vs_alphaL), max(sw.sep_vs_alphaL))
 'Vulnerability exponent', ...
   sprintf('%g to %g', min(sw.gammaV), max(sw.gammaV)), '-', ...
   sprintf('median separation %.0f to %.0f m, goal advantage %.2f to %.2f%%', ...
           min(sw.sep_vs_gamma), max(sw.sep_vs_gamma), ...
           100*min(sw.adv_vs_gamma), 100*max(sw.adv_vs_gamma))
 };
t = local_mdtable({'Assumption varied', 'Range', 'R', 'Other effect'}, rows);
end

% =========================================================================
function t = local_mdtable(head, rows)
lines = {['| ' strjoin(head, ' | ') ' |'], ...
         ['|' repmat('---|', 1, numel(head))]};
for k = 1:size(rows, 1)
    lines{end+1} = ['| ' strjoin(rows(k, :), ' | ') ' |']; %#ok<AGROW>
end
t = strjoin(lines, newline);
end
