function fk_write_results(cfg, R)
%FK_WRITE_RESULTS  Write results/results.json, the only source the manuscript
%                  quotes numbers from.

r = struct();
r.run_date = R.run_date;
r.matlab_version = R.matlab_version;
r.runtime_s = R.runtime_s;
r.seed = cfg.seed;

% ------------------------------------------------------------------ data ---
r.data.n_wells_total      = R.D.n_total;
r.data.n_wells_inside     = R.D.n_inside;
r.data.n_with_K           = R.D.n_K;
r.data.n_with_NO3         = R.D.n_NO3;
r.data.n_inside_with_K    = R.D.n_inside_K;
r.data.n_inside_with_NO3  = R.D.n_inside_NO3;
r.data.n_boundary_vertices = R.D.n_boundary;
r.data.K_min  = min(R.D.K(R.D.has_K));
r.data.K_max  = max(R.D.K(R.D.has_K));
r.data.K_mean = mean(R.D.K(R.D.has_K));
r.data.K_sd   = std(R.D.K(R.D.has_K));
r.data.NO3_min  = min(R.D.NO3(R.D.has_NO3));
r.data.NO3_max  = max(R.D.NO3(R.D.has_NO3));
r.data.NO3_mean = mean(R.D.NO3(R.D.has_NO3));
r.data.NO3_sd   = std(R.D.NO3(R.D.has_NO3));
r.data.DRASTIC_min  = min(R.D.DRASTIC);
r.data.DRASTIC_max  = max(R.D.DRASTIC);
r.data.DRASTIC_mean = mean(R.D.DRASTIC);
r.data.DRASTIC_sd   = std(R.D.DRASTIC);
r.data.facies_counts = [sum(R.D.has_K & R.D.facies == 1), ...
                        sum(R.D.has_K & R.D.facies == 2), ...
                        sum(R.D.has_K & R.D.facies == 3)];
r.data.facies_K_eta2 = R.P.eta2_facies_K;
r.data.facies_K_F    = R.P.F_facies_K;
r.data.facies_K_p    = R.P.p_facies_K;

% -------------------------------------------------------- predictability ---
r.predictability.n = R.Pr.n;
r.predictability.idw_power = R.Pr.idw_power;
r.predictability.idw_R2 = R.Pr.idw_R2;
r.predictability.kriging_range_m = R.Pr.krig_range;
r.predictability.kriging_R2 = R.Pr.krig_R2;
r.predictability.linear_names = R.Pr.lin_name;
r.predictability.linear_R2 = R.Pr.lin_R2;
r.predictability.best_R2 = R.Pr.best_R2;
r.predictability.corr_northing = R.Pr.corr_northing;
r.predictability.corr_easting = R.Pr.corr_easting;
r.predictability.corr_alongflow = R.Pr.corr_alongflow;
r.predictability.corr_drastic = R.Pr.corr_drastic;
r.predictability.corr_K = R.Pr.corr_K;
r.predictability.corr_depth = R.Pr.corr_depth;

% ---------------------------------------------------------- grid, prior ----
r.grid.dx_m = cfg.dx;
r.grid.nx = R.G.nx;  r.grid.ny = R.G.ny;
r.grid.n_active = R.G.n_active;
r.grid.area_km2 = R.G.area_km2;
r.grid.area_polygon_km2 = R.G.area_polygon_km2;
r.grid.flow_azimuth_deg = cfg.azim_flow_deg;

r.prior.sd_independent = R.P.sd_iv;
r.prior.sd_facies_correlated = R.P.sd_facies;
r.prior.rho_used = R.P.rho_used;
r.prior.sd_used = R.P.sd_facies(2);
r.prior.sd_inflation_pct = 100 * (R.P.sd_facies(2) / R.P.sd_iv - 1);
r.prior.network_mean = R.P.ne_network_mean;
r.prior.network_sd = R.P.ne_network_sd;
r.prior.ne_grid_min = min(R.P.ne);
r.prior.ne_grid_max = max(R.P.ne);
r.prior.ne_grid_mean = mean(R.P.ne);
r.prior.variogram_nugget = R.P.nugget;
r.prior.variogram_sill = R.P.sill;
r.prior.variogram_range_m = R.P.vrange;
r.prior.n_kl = numel(R.P.lam);
r.prior.kl_variance_captured = R.P.kl_var_captured;

% --------------------------------------------------------------- forward ---
r.forward.gradient = R.F.i_gradient;
r.forward.darcy_flux_m_per_day = R.F.q0;
r.forward.alphaL_m = cfg.alphaL;
r.forward.alphaT_m = cfg.alphaT_ratio * cfg.alphaL;
r.forward.Gamma_min = min(R.F.Gam);
r.forward.Gamma_max = max(R.F.Gam);

% ----------------------------------------------------------- calibration ---
r.calibration.n_obs = R.Cal.n_obs;
r.calibration.basis = R.Cal.psi_names;
r.calibration.T0_days = R.Cal.T0;
r.calibration.T0_years = R.Cal.T0 / 365.25;
r.calibration.rmse_mg_per_L = R.Cal.rmse;
r.calibration.R2 = R.Cal.R2;
r.calibration.R2_without_DRASTIC = R.Cal.R2_noV;
r.calibration.loo_rmse_mg_per_L = R.Cal.loo_rmse;
r.calibration.loo_R2 = R.Cal.loo_R2;
r.calibration.modelled_min = min(R.Cal.C);
r.calibration.modelled_max = max(R.Cal.C);
r.calibration.frac_source_negative_before_clipping = R.Cal.frac_negative;
r.calibration.null_kappa = R.Cal.null_kappa;
r.calibration.null_R2 = R.Cal.null_R2;
r.calibration.null_best_R2 = R.Cal.null_best_R2;
r.calibration.null_best_kappa = R.Cal.null_best_kappa;
r.calibration.null_best_T0_days = R.Cal.null_best_T0;

% --------------------------------------------------- Fisher information ----
r.fisher.trace_prior = R.S.tr_prior;
r.fisher.trace_observations = R.S.tr_obs;
r.fisher.observation_share = R.S.frac_obs;
r.fisher.dofs = R.S.dofs;
r.fisher.n_modes = R.S.n_modes;
r.fisher.logdet_prior = R.S.logdet_prior;
r.fisher.logdet_total = R.S.logdet_tot;
r.fisher.information_gain_nats = R.S.info_gain_nats;
r.fisher.sd_reduction_fraction = R.S.sd_reduction;
r.fisher.finite_difference_check = R.S.fd_check_rel;
r.fisher.fd_steps = R.S.fd_steps;
r.fisher.fd_relative = R.S.fd_rel;
r.fisher.fd_convergence_order = R.S.fd_slope;
r.fisher.sd_prior_mean = mean(R.S.sd_prior_field);
r.fisher.sd_post_mean = mean(R.S.sd_post);

% ---------------------------------------------------------------- design ---
r.design.n_candidates = R.Des.n_candidates;
r.design.n_new_wells = cfg.n_new_wells;
r.design.sd_local = cfg.sd_local;
r.design.gamma_V = R.Des.gammaV;
r.design.tau_weighted_mean_days = R.Des.tau_mean_prior;
r.design.sd_goal_prior_days = R.Des.sd_goal_prior;
r.design.sd_goal_after_D_days = sqrt(R.Des.D.var_goal(end));
r.design.sd_goal_after_goal_days = sqrt(R.Des.G.var_goal(end));
r.design.D_cells_xy = [R.G.xa(R.Des.D.cells), R.G.ya(R.Des.D.cells)];
r.design.goal_cells_xy = [R.G.xa(R.Des.G.cells), R.G.ya(R.Des.G.cells)];
r.design.D_DRASTIC = R.Des.V_D;
r.design.goal_DRASTIC = R.Des.V_G;
r.design.domain_mean_DRASTIC = mean(R.P.V);
r.design.separation_m = R.Des.separation_m;
r.design.median_separation_m = median(R.Des.separation_m);
r.design.logdet_gain_D = R.Des.D.logdet_gain;
r.design.information_horizon_m = R.Des.horizon;
r.design.information_correlation_max = R.Des.info_corr_max;
r.design.information_correlation_unconstrained_max = R.Des.info_corr_free_max;
r.design.information_correlation_horizon_max = R.Des.info_corr_hor_max;
r.design.information_horizon_median_m = R.Des.horizon_median;
r.design.min_separation_horizon_m = R.Des.min_sep_horizon;
r.design.sd_goal_after_D_horizon_days = sqrt(R.Des.D_hor.var_goal(end));
r.design.sd_goal_after_goal_horizon_days = sqrt(R.Des.G_hor.var_goal(end));
r.design.unconstrained_span_m = R.Des.free_span_m;
r.design.min_separation_m = R.Des.min_sep;
r.design.horizon_fraction = R.Des.horizon_frac;
r.design.gain_conc_best = R.Des.gain_conc_best;
r.design.gain_pump_best = R.Des.gain_pump_best;
r.design.sd_goal_after_D_free_days = sqrt(R.Des.D_free.var_goal(end));
r.design.sd_goal_after_goal_free_days = sqrt(R.Des.G_free.var_goal(end));
r.design.R_equivalence = R.Des.R_equiv;
r.design.sd_pump_test_on_ne = R.Des.sd_pump_ne;

% --------------------------------------------------------- local operator --
r.local.var_reduction_D = R.Loc.var_reduction_D;
r.local.var_reduction_goal = R.Loc.var_reduction_goal;
r.local.goal_advantage = R.Loc.goal_advantage;
r.local.median_separation_m = R.Loc.median_separation_m;

% ------------------------------------------------------------- sampling ----
r.sampling.gn_iterations = R.Par.gn_iters_used;
r.sampling.F_at_map = R.Par.phi_F;
r.sampling.condition_number = R.Par.cond_H;
r.sampling.iters_natural = R.Par.iters_natural;
r.sampling.iters_euclidean = R.Par.iters_euclidean;
r.sampling.speedup = R.Par.speedup;
r.sampling.sqrt_condition = sqrt(R.Par.cond_H);
r.sampling.M_implicit = cfg.M_particles;
r.sampling.ess_implicit = R.Par.ess_implicit;
r.sampling.M_sir = cfg.M_sir;
r.sampling.ess_sir = R.Par.ess_sir;
r.sampling.n_forward_solves = R.Par.n_forward;

% ------------------------------------------------------------------ MRVA ---
r.mrva.n_wells = R.M.n_wells;
r.mrva.drastic_mismatch = R.M.drastic_mismatch;
r.mrva.K_min = R.M.K_range(1);  r.mrva.K_max = R.M.K_range(2);
r.mrva.K_mean = R.M.K_mean;     r.mrva.K_sd = R.M.K_sd;
r.mrva.V_min = R.M.V_range(1);  r.mrva.V_max = R.M.V_range(2);
r.mrva.V_mean = R.M.V_mean;     r.mrva.V_sd = R.M.V_sd;
r.mrva.variogram_range_m = R.M.range;
r.mrva.area_km2 = R.M.area_km2;
r.mrva.n_active = R.M.n_active;
r.mrva.var_reduction_D = R.M.oed.var_reduction_D;
r.mrva.var_reduction_goal = R.M.oed.var_reduction_goal;
r.mrva.goal_advantage = R.M.oed.goal_advantage;
r.mrva.median_separation_m = R.M.oed.median_separation_m;
r.mrva.mean_DRASTIC_D = mean(R.M.V_D);
r.mrva.mean_DRASTIC_goal = mean(R.M.V_goal);

% ----------------------------------------------------------------- sweep ---
r.sweep = R.Sw;

txt = jsonencode(r, 'PrettyPrint', true);
fid = fopen(fullfile(cfg.dir_results, 'results.json'), 'w');
fwrite(fid, txt);
fclose(fid);
fprintf('[out] results/results.json written (%d characters)\n', numel(txt));
end
