function fk_figures(R)
%FK_FIGURES  Every figure of the manuscript, from results/run_all.mat.
%
%  Figures are drawn at the width they are printed at (90 mm one column,
%  190 mm two columns), with Position equal to PaperPosition so that nothing
%  is rescaled on export, and written at 600 dpi as PNG and TIFF.

if nargin < 1
    cfg0 = fk_config();
    R = load(fullfile(cfg0.dir_results, 'run_all.mat'));
end
cfg = R.cfg;  D = R.D;  G = R.G;  P = R.P;  Cal = R.Cal;  S = R.S;
Des = R.Des;  Loc = R.Loc;  Par = R.Par;  M = R.M;  Sw = R.Sw;  Pr = R.Pr;
if ~exist(cfg.dir_figures, 'dir'), mkdir(cfg.dir_figures); end

set(0, 'DefaultAxesFontName', 'Arial', 'DefaultTextFontName', 'Arial', ...
       'DefaultAxesFontSize', 8, 'DefaultTextFontSize', 8, ...
       'DefaultAxesLineWidth', 0.5, 'DefaultLineLineWidth', 0.9);
W2 = 19.0;
CB = [0.00 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19; 0.49 0.18 0.56];

% ============================================================== Figure 1 ===
f = newfig(W2, 10.0);
subplot(1, 3, 1);
plot(D.bx/1000, D.by/1000, 'k-', 'LineWidth', 0.6); hold on
plot(D.x(~D.has_K)/1000, D.y(~D.has_K)/1000, '.', 'Color', [.65 .65 .65], ...
     'MarkerSize', 4);
plot(D.x(D.has_K)/1000, D.y(D.has_K)/1000, 'o', 'MarkerSize', 3, ...
     'MarkerFaceColor', CB(1,:), 'MarkerEdgeColor', 'none');
plot(D.x(D.has_NO3)/1000, D.y(D.has_NO3)/1000, '^', 'MarkerSize', 3.4, ...
     'MarkerFaceColor', CB(2,:), 'MarkerEdgeColor', 'none');
axis equal tight; box on
xlabel('Easting (km)'); ylabel('Northing (km)');
title(sprintf('(a) %d points, %d with K, %d with NO_3', ...
              D.n_total, D.n_K, D.n_NO3), 'FontWeight', 'normal');
legend({'boundary', 'DRASTIC only', 'pumping test', 'nitrate'}, ...
       'Location', 'southoutside', 'Box', 'off');

subplot(1, 3, 2);
scatterfield(G, P.V); hold on
plot(D.bx/1000, D.by/1000, 'k-', 'LineWidth', 0.4);
quiverflow(G);
c = colorbar('southoutside'); c.Label.String = 'DRASTIC index';
title('(b) vulnerability and regional flow', 'FontWeight', 'normal');

subplot(1, 3, 3);
scatterfield(G, P.ne);
c = colorbar('southoutside'); c.Label.String = 'effective porosity prior';
title('(c) porosity prior', 'FontWeight', 'normal');
printfig(f, cfg, 'fig1_study_area');

% ============================================================== Figure 2 ===
f = newfig(W2, 7.5);
subplot(1, 2, 1);
names = [arrayfun(@(p) sprintf('IDW %d', p), Pr.idw_power, 'uni', 0), ...
         arrayfun(@(r) sprintf('kriging %g km', r/1000), Pr.krig_range, 'uni', 0), ...
         Pr.lin_name, {'transport model'}];
vals = [Pr.idw_R2, Pr.krig_R2, Pr.lin_R2, Cal.loo_R2];
barh(vals, 'FaceColor', CB(1,:), 'EdgeColor', 'none'); hold on
set(gca, 'YTick', 1:numel(vals), 'YTickLabel', names, 'YDir', 'reverse', ...
    'TickLabelInterpreter', 'none');
xline(0, 'k-', 'LineWidth', 0.9);
xlabel('leave-one-out R^2'); box on
title('(a) nothing predicts a held-out well', 'FontWeight', 'normal');

subplot(1, 2, 2);
plot(Cal.y_obs, Cal.y_fit, 'o', 'MarkerSize', 4, ...
     'MarkerFaceColor', CB(1,:), 'MarkerEdgeColor', 'none'); hold on
plot(Cal.y_obs, Cal.loo_pred, 's', 'MarkerSize', 4, ...
     'MarkerFaceColor', CB(2,:), 'MarkerEdgeColor', 'none');
lim = [0 max([Cal.y_obs; Cal.y_fit; Cal.loo_pred]) * 1.15];
plot(lim, lim, 'k--', 'LineWidth', 0.6);
xlim(lim); ylim(lim); axis square; box on
xlabel('observed NO_3 (mg L^{-1})'); ylabel('modelled NO_3 (mg L^{-1})');
legend({sprintf('in sample, R^2 = %.2f', Cal.R2), ...
        sprintf('leave one out, R^2 = %.2f', Cal.loo_R2), '1:1'}, ...
       'Location', 'northwest', 'Box', 'off');
title('(b) in sample against out of sample', 'FontWeight', 'normal');
printfig(f, cfg, 'fig2_predictability');

% ============================================================== Figure 3 ===
f = newfig(W2, 7.0);
subplot(1, 3, 1);
semilogy(P.lam / P.lam(1), 'k-'); box on
xlabel('Karhunen-Loeve mode'); ylabel('\lambda_k / \lambda_1');
title(sprintf('(a) %d modes, %.0f%% of prior variance', ...
              numel(P.lam), 100*P.kl_var_captured), 'FontWeight', 'normal');
subplot(1, 3, 2);
plot(P.vario_h/1000, P.vario_g, 'o', 'MarkerSize', 3.5, ...
     'MarkerFaceColor', CB(1,:), 'MarkerEdgeColor', 'none'); hold on
hh = linspace(0, max(P.vario_h), 200);
plot(hh/1000, P.nugget + P.sill*(1 - exp(-hh/P.vrange)), 'k-');
xlabel('lag (km)'); ylabel('\gamma of log_{10}K'); box on
title(sprintf('(b) exponential range %.0f m', P.vrange), 'FontWeight', 'normal');
subplot(1, 3, 3);
scatterfield(G, Cal.source / (cfg.dx*cfg.dy));
c = colorbar('southoutside'); c.Label.String = 'loading (mg L^{-1} d^{-1})';
title(sprintf('(c) fitted loading, T_0 = %.0f yr', Cal.T0/365.25), ...
      'FontWeight', 'normal');
printfig(f, cfg, 'fig3_prior_and_loading');

% ============================================================== Figure 4 ===
f = newfig(W2, 7.0);
subplot(1, 3, 1);
scatterfield(G, Cal.C);
c = colorbar('southoutside'); c.Label.String = 'NO_3 (mg L^{-1})';
title('(a) modelled nitrate at T_0', 'FontWeight', 'normal');
subplot(1, 3, 2);
scatterfield(G, S.sd_prior_field);
c = colorbar('southoutside'); c.Label.String = 'prior sd of n_e';
title('(b) prior standard deviation', 'FontWeight', 'normal');
subplot(1, 3, 3);
scatterfield(G, S.sd_prior_field - S.sd_post);
c = colorbar('southoutside'); c.Label.String = 'sd reduction';
title(sprintf('(c) what the %d nitrate wells buy', Cal.n_obs), ...
      'FontWeight', 'normal');
printfig(f, cfg, 'fig4_information');

% ============================================================== Figure 5 ===
% The transport-operator design, and the greedy clustering pathology.
f = newfig(W2, 9.0);
subplot(1, 3, 1);
scatterfield(G, log10(max(Des.gain_conc, realmin))); hold on
plot(G.xa(Des.D_free.cells)/1000, G.ya(Des.D_free.cells)/1000, 'ko', ...
     'MarkerSize', 5, 'LineWidth', 1);
c = colorbar('southoutside'); c.Label.String = 'log_{10} variance gain';
title(sprintf('(a) unconstrained greedy: %d wells in %.0f m', ...
              cfg.n_new_wells, Des.free_span_m), 'FontWeight', 'normal');

subplot(1, 3, 2);
scatterfield(G, P.V); hold on
plot(G.xa(Des.D.cells)/1000, G.ya(Des.D.cells)/1000, 'ks', ...
     'MarkerSize', 6, 'LineWidth', 1);
plot(G.xa(Des.G.cells)/1000, G.ya(Des.G.cells)/1000, 'rp', ...
     'MarkerSize', 8, 'LineWidth', 1);
c = colorbar('southoutside'); c.Label.String = 'DRASTIC index';
title(sprintf('(b) separated by %.0f m', Des.min_sep), 'FontWeight', 'normal');

subplot(1, 3, 3);
plot(0:cfg.n_new_wells, sqrt(Des.D.var_goal), 'ks-', 'MarkerSize', 4, ...
     'MarkerFaceColor', 'k'); hold on
plot(0:cfg.n_new_wells, sqrt(Des.G.var_goal), 'r^--', 'MarkerSize', 4, ...
     'MarkerFaceColor', 'r');
xlabel('concentration wells added');
ylabel('sd of weighted arrival time (d)');
legend({'D-optimal', 'goal-oriented'}, 'Box', 'off', 'Location', 'northeast');
box on; title('(c) almost nothing to allocate', 'FontWeight', 'normal');
printfig(f, cfg, 'fig5_transport_design');

% ============================================================== Figure 6 ===
% The local-operator design, where the criterion does matter, at both sites.
f = newfig(W2, 9.0);
subplot(1, 3, 1);
scatterfield(G, P.V); hold on
plot(G.xa(Loc.D.cells)/1000, G.ya(Loc.D.cells)/1000, 'ks', ...
     'MarkerSize', 6, 'LineWidth', 1);
plot(G.xa(Loc.goal.cells)/1000, G.ya(Loc.goal.cells)/1000, 'rp', ...
     'MarkerSize', 8, 'LineWidth', 1);
c = colorbar('southoutside'); c.Label.String = 'DRASTIC index';
title('(a) Fushe-Kuqe, characterisation points', 'FontWeight', 'normal');

subplot(1, 3, 2);
nw = cfg.n_new_wells;
plot(0:nw, Loc.D.var_goal / Loc.var_prior, 'ks-', 'MarkerSize', 4, ...
     'MarkerFaceColor', 'k'); hold on
plot(0:nw, Loc.goal.var_goal / Loc.var_prior, 'r^--', 'MarkerSize', 4, ...
     'MarkerFaceColor', 'r');
plot(0:nw, M.oed.D.var_goal / M.oed.var_prior, 'ks:', 'MarkerSize', 4);
plot(0:nw, M.oed.goal.var_goal / M.oed.var_prior, 'r^:', 'MarkerSize', 4);
xlabel('points added'); ylabel('Var(g) relative to prior'); ylim([0 1]);
legend({'Fushe-Kuqe D', 'Fushe-Kuqe goal', 'MRVA D', 'MRVA goal'}, ...
       'Box', 'off', 'Location', 'southwest');
box on; title('(b) goal orientation pays', 'FontWeight', 'normal');

subplot(1, 3, 3);
scatter(M.xa/1000, M.ya/1000, 4, M.V, 'filled', 'Marker', 's'); hold on
plot(M.xa(M.oed.D.cells)/1000, M.ya(M.oed.D.cells)/1000, 'ks', ...
     'MarkerSize', 6, 'LineWidth', 1);
plot(M.xa(M.oed.goal.cells)/1000, M.ya(M.oed.goal.cells)/1000, 'rp', ...
     'MarkerSize', 8, 'LineWidth', 1);
axis equal tight; box on
xlabel('x (km)'); ylabel('y (km)');
c = colorbar('southoutside'); c.Label.String = 'DRASTIC index';
title(sprintf('(c) MRVA, %d USGS wells', M.n_wells), 'FontWeight', 'normal');
printfig(f, cfg, 'fig6_local_design_and_replication');

% ============================================================== Figure 7 ===
f = newfig(W2, 6.5);
subplot(1, 4, 1);
bar([Par.ess_implicit, Par.ess_sir], 'FaceColor', CB(1,:), 'EdgeColor', 'none');
set(gca, 'XTickLabel', {sprintf('implicit\nM=%d', cfg.M_particles), ...
                        sprintf('SIR\nM=%d', cfg.M_sir)});
ylabel('normalised ESS'); box on; ylim([0 1]);
title('(a) both samplers are easy', 'FontWeight', 'normal');

subplot(1, 4, 2);
plot(Sw.sigma_obs, Sw.R_vs_sigma, 'ko-', 'MarkerFaceColor', 'k', ...
     'MarkerSize', 4); hold on
yline(1, 'r--');
xlabel('\sigma_{obs} (mg L^{-1})'); ylabel('R'); box on
title('(b) R against noise', 'FontWeight', 'normal');

subplot(1, 4, 3);
semilogx(Sw.alphaL, Sw.R_vs_alphaL, 'ko-', 'MarkerFaceColor', 'k', ...
         'MarkerSize', 4); hold on
yline(1, 'r--');
xlabel('\alpha_L (m)'); ylabel('R'); box on
title('(c) R against dispersivity', 'FontWeight', 'normal');

subplot(1, 4, 4);
yyaxis left;  plot(Sw.gammaV, Sw.sep_vs_gamma/1000, 'o-', 'MarkerSize', 4);
ylabel('median separation (km)');
yyaxis right; plot(Sw.gammaV, 100*Sw.adv_vs_gamma, 's-', 'MarkerSize', 4);
ylabel('goal advantage (%)'); xlabel('\gamma_V'); box on
title('(d) weighting exponent', 'FontWeight', 'normal');
printfig(f, cfg, 'fig7_sampling_and_sensitivity');

% ==================================================== graphical abstract ===
% Elsevier asks for 531 x 1328 px at 96 dpi minimum; 13.3 x 5.3 cm at 600 dpi
% is comfortably above that.
f = newfig(13.3, 5.3);
ax1 = axes('Position', [0.04 0.18 0.28 0.66]);
scatterfield(G, P.V); hold on
plot(D.x(D.has_NO3)/1000, D.y(D.has_NO3)/1000, '^', 'MarkerSize', 3.2, ...
     'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'none');
set(ax1, 'XTickLabel', [], 'YTickLabel', []);
xlabel(''); ylabel('');
title(sprintf('%d nitrate wells', Cal.n_obs), 'FontWeight', 'normal', ...
      'FontSize', 8);

ax2 = axes('Position', [0.345 0.12 0.265 0.76]);
% Typography rather than a chart: the quantity is a fraction of 150 so small
% that no shared axis can show it and be read at abstract size.
axis off; xlim([0 1]); ylim([0 1]); hold on
text(0.5, 0.92, 'the network already has', 'HorizontalAlignment', 'center', ...
     'FontSize', 8.5);
text(0.5, 0.74, sprintf('%.3f', S.dofs), 'HorizontalAlignment', 'center', ...
     'FontSize', 26, 'FontWeight', 'bold', 'Color', CB(2,:));
text(0.5, 0.57, sprintf('of %d directions constrained', numel(P.lam)), ...
     'HorizontalAlignment', 'center', 'FontSize', 8.5);
plot([0.12 0.88], [0.46 0.46], 'k-', 'LineWidth', 0.6);
text(0.5, 0.35, 'one nitrate sample buys', 'HorizontalAlignment', 'center', ...
     'FontSize', 8.5);
text(0.5, 0.20, sprintf('%.2f', Des.R_equiv), 'HorizontalAlignment', 'center', ...
     'FontSize', 22, 'FontWeight', 'bold', 'Color', CB(2,:));
text(0.5, 0.05, 'of a pumping test', 'HorizontalAlignment', 'center', ...
     'FontSize', 8.5);

ax3 = axes('Position', [0.735 0.24 0.245 0.58]);
plot(0:cfg.n_new_wells, Loc.D.var_goal / Loc.var_prior, 'ks-', ...
     'MarkerSize', 4, 'MarkerFaceColor', 'k'); hold on
plot(0:cfg.n_new_wells, Loc.goal.var_goal / Loc.var_prior, 'r^--', ...
     'MarkerSize', 4, 'MarkerFaceColor', 'r');
xlabel('characterisation points'); ylabel('Var(g) / prior');
ylim([0 1]); box on
legend({'D-optimal', 'goal-oriented'}, 'Box', 'off', 'FontSize', 7, ...
       'Location', 'southwest');
title('where design pays', 'FontWeight', 'normal', 'FontSize', 8);
printfig(f, cfg, 'graphical_abstract');

fprintf('[fig] seven figures and a graphical abstract written to %s\n', ...
        cfg.dir_figures);
end

% =========================================================================
function f = newfig(wcm, hcm)
f = figure('Units', 'centimeters', 'Position', [2 2 wcm hcm], 'Color', 'w');
set(f, 'PaperUnits', 'centimeters', 'PaperSize', [wcm hcm], ...
       'PaperPosition', [0 0 wcm hcm], 'PaperPositionMode', 'manual');
end

function printfig(f, cfg, name)
print(f, fullfile(cfg.dir_figures, [name '.png']), '-dpng', '-r600');
print(f, fullfile(cfg.dir_figures, [name '.tif']), '-dtiff', '-r600');
close(f);
end

function scatterfield(G, v)
scatter(G.xa/1000, G.ya/1000, 1.5, v, 'filled', 'Marker', 's');
axis equal tight; box on
xlabel('Easting (km)'); ylabel('Northing (km)');
end

function quiverflow(G)
xc = mean(G.xa)/1000;  yc = mean(G.ya)/1000;
L = 0.15 * (max(G.xa) - min(G.xa)) / 1000;
quiver(xc, yc, G.exi(1)*L, G.exi(2)*L, 0, 'k', 'LineWidth', 1.2, ...
       'MaxHeadSize', 2);
text(xc + G.exi(1)*L, yc + G.exi(2)*L, ' flow', 'FontSize', 7);
end
