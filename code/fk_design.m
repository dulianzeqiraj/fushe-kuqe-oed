function Des = fk_design(cfg, D, G, P, F, Cal, S)
%FK_DESIGN  Classical D-optimal and goal-oriented monitoring network design.
%
%  Two criteria are compared.
%
%  (1) Classical D-optimality.  Adding one concentration observation at x
%      updates the Fisher information by the rank-one term sigma^-2 J(x)'J(x),
%      and by the matrix determinant lemma
%          log det I+ - log det I = log(1 + sigma^-2 J(x) Sigma J(x)').
%      This is the objective maximised by the classical criterion, and it is
%      monotone submodular in the chosen set, so the greedy sequence carries
%      the standard (1 - 1/e) guarantee (Krause et al. 2008).
%
%  (2) Goal-oriented design for the vulnerability-weighted arrival time.
%      Protection-zone delineation does not depend on the porosity field as
%      such; it depends on how long a contaminant takes to reach the discharge
%      boundary, and it depends on that most where the aquifer is vulnerable.
%      With a uniform Darcy flux the residual travel time from x is
%      tau(x) = L(x) n_e(x) / q0, where L(x) is the remaining distance along
%      the flow direction, so a porosity error maps linearly onto an arrival
%      time error.  Define the vulnerability-weighted mean arrival-time error
%          g(theta) = sum_x w(x) L(x) dn_e(x) / (q0 sum_x w(x)),
%          w(x) = (V(x)/mean V)^gamma,
%      a linear functional g = c'theta.  The design minimises Var(g) = c'Sigma c,
%      and one observation at x reduces it by
%          Delta_g(x) = sigma^-2 (c' Sigma J(x)')^2 / (1 + sigma^-2 J(x) Sigma J(x)').
%
%  Note on the guarantee.  Submodularity is a property of the log-determinant
%  objective.  It is not claimed for the goal-oriented variance reduction,
%  where the greedy sequence is used as a heuristic and reported as such.

na = G.n_active;
sig2 = cfg.sigma_obs^2;

% ------------------------------------- candidate set: every active cell ----
% Sites already occupied by a nitrate well are excluded.
cand = true(na, 1);
cand(Cal.cell_of) = false;
Des.n_candidates = nnz(cand);

% ------------------------------------------- goal functional coefficients --
% Remaining distance to the downstream boundary along the flow direction.
xi_a = G.xa * G.exi(1) + G.ya * G.exi(2);
Lrem = max(xi_a) - xi_a;
w = (P.V / mean(P.V)) .^ cfg.gamma_Vdesign;
cw = w .* Lrem;
c_goal = (P.Phi' * cw) / (F.q0 * sum(w));           % days per unit theta
Des.c_goal = c_goal;
Des.Lrem = Lrem;
Des.w = w;
Des.tau_mean_prior = sum(w .* Lrem .* P.ne) / (F.q0 * sum(w));
fprintf(['[design] vulnerability-weighted mean residual travel time ' ...
         '%.0f d (%.1f yr)\n'], Des.tau_mean_prior, Des.tau_mean_prior/365.25);

% ------------------------------------------------- information horizon -----
% Computed before the placement, because it is what the minimum separation
% should be set from.  The kernel through which one observation sees the field
% is kappa_x(x') = sum_k lambda_k J(x,k) phi_k(x'), the prior covariance
% between the concentration at x and the porosity elsewhere; the horizon is the
% radius holding a given fraction of its absolute mass.  Reference sites are
% the single best one-shot gain and the first three existing nitrate wells.
one_shot = sum((S.J * S.Sigma) .* S.J, 2);
one_shot(~cand) = -inf;
[~, best_single] = max(one_shot);
Des.horizon_sites = [best_single; Cal.cell_of(1:3)];
Des.horizon = zeros(numel(Des.horizon_sites), 1);
for t = 1:numel(Des.horizon_sites)
    ci = Des.horizon_sites(t);
    a = abs(P.Phi * (P.lam .* S.J(ci, :)'));
    rr = hypot(G.xa - G.xa(ci), G.ya - G.ya(ci));
    [rs, ord] = sort(rr);
    cum = cumsum(a(ord)) / sum(a);
    Des.horizon(t) = rs(find(cum >= cfg.horizon_frac, 1));
end
Des.horizon_median = median(Des.horizon);
fprintf(['[design] information horizon at %.0f%% of kernel mass: %.0f - %.0f m ' ...
         '(median %.0f m), against a variogram range of %.0f m\n'], ...
        100*cfg.horizon_frac, min(Des.horizon), max(Des.horizon), ...
        Des.horizon_median, P.vrange);

% ------------------------------------------------------ greedy placement ---
% Unconstrained greedy first, then the same search under two minimum
% separations.  The unconstrained version is kept because what it does is a
% result in its own right: on a 100 m grid whose information field varies over
% kilometres, the rank-one update removes so little of the information at the
% chosen cell that the next choice is its neighbour, and the sequence collapses
% into a cluster of mutually redundant wells.  The two constraints test which
% length scale has to be forbidden, the variogram range or the information
% horizon.
Des.gammaV = cfg.gamma_Vdesign;
Des.min_sep = P.vrange;
Des.min_sep_horizon = Des.horizon_median;
Des.D_hor = local_greedy(S.Sigma, S.J, c_goal, cand, sig2, cfg.n_new_wells, ...
                         'D', Des.min_sep_horizon, G);
Des.G_hor = local_greedy(S.Sigma, S.J, c_goal, cand, sig2, cfg.n_new_wells, ...
                         'goal', Des.min_sep_horizon, G);
Des.D_free = local_greedy(S.Sigma, S.J, c_goal, cand, sig2, cfg.n_new_wells, ...
                          'D', 0, G);
Des.G_free = local_greedy(S.Sigma, S.J, c_goal, cand, sig2, cfg.n_new_wells, ...
                          'goal', 0, G);
Des.D  = local_greedy(S.Sigma, S.J, c_goal, cand, sig2, cfg.n_new_wells, ...
                      'D', Des.min_sep, G);
Des.G  = local_greedy(S.Sigma, S.J, c_goal, cand, sig2, cfg.n_new_wells, ...
                      'goal', Des.min_sep, G);

Des.free_span_m = max(hypot(G.xa(Des.D_free.cells) - G.xa(Des.D_free.cells)', ...
                            G.ya(Des.D_free.cells) - G.ya(Des.D_free.cells)'), ...
                      [], 'all');
fprintf(['[design] unconstrained greedy D-optimal spreads its %d wells over ' ...
         '%.0f m; minimum separation set to %.0f m\n'], ...
        cfg.n_new_wells, Des.free_span_m, Des.min_sep);

for nm = {'D', 'G'}
    r = Des.(nm{1});
    fprintf('[design] %s-optimal sequence:\n', nm{1});
    for j = 1:numel(r.cells)
        ci = r.cells(j);
        fprintf(['   %d: x %.0f y %.0f  DRASTIC %.0f  ' ...
                 'logdet gain %.4f  Var(g) %.4g -> sd(g) %.1f d\n'], ...
                j, G.xa(ci), G.ya(ci), P.V(ci), r.logdet_gain(j), ...
                r.var_goal(j + 1), sqrt(r.var_goal(j + 1)));
    end
end

Des.var_goal_prior = c_goal' * S.Sigma * c_goal;
Des.sd_goal_prior = sqrt(Des.var_goal_prior);
fprintf(['[design] sd of the vulnerability-weighted arrival time: ' ...
         '%.1f d before augmentation\n'], Des.sd_goal_prior);
fprintf(['[design] after 5 wells: D-optimal %.1f d, goal-oriented %.1f d ' ...
         '(%.1f%% better)\n'], sqrt(Des.D.var_goal(end)), ...
        sqrt(Des.G.var_goal(end)), ...
        100 * (1 - sqrt(Des.G.var_goal(end)) / sqrt(Des.D.var_goal(end))));

% separation between the two designs
dd = zeros(cfg.n_new_wells, 1);
for j = 1:cfg.n_new_wells
    dd(j) = min(hypot(G.xa(Des.G.cells(j)) - G.xa(Des.D.cells), ...
                      G.ya(Des.G.cells(j)) - G.ya(Des.D.cells)));
end
Des.separation_m = dd;
Des.V_D = P.V(Des.D.cells);
Des.V_G = P.V(Des.G.cells);
fprintf(['[design] mean DRASTIC at the chosen sites: D-optimal %.1f, ' ...
         'goal-oriented %.1f (domain mean %.1f)\n'], ...
        mean(Des.V_D), mean(Des.V_G), mean(P.V));
fprintf('[design] median distance between the two designs %.0f m\n', median(dd));

% --------------------------------------- information exclusion diagnostic --
% Correlation of the information carried by two observations.  A greedy
% sequence that keeps this small is placing non-redundant wells.
nn = cfg.n_new_wells;
    function Rm = infocorr(cells)
        Rm = eye(nn);
        for a = 1:nn
            for b = a+1:nn
                Ja = S.J(cells(a), :);  Jb = S.J(cells(b), :);
                Rm(a, b) = (Ja * S.Sigma * Jb') / ...
                           sqrt((Ja * S.Sigma * Ja') * (Jb * S.Sigma * Jb'));
                Rm(b, a) = Rm(a, b);
            end
        end
    end
Des.info_corr = infocorr(Des.G.cells);
Des.info_corr_free = infocorr(Des.G_free.cells);
Des.info_corr_hor = infocorr(Des.G_hor.cells);
mx = @(Rm) max(abs(Rm(~eye(nn))));
Des.info_corr_max = mx(Des.info_corr);
Des.info_corr_free_max = mx(Des.info_corr_free);
Des.info_corr_hor_max = mx(Des.info_corr_hor);
fprintf(['[design] maximum pairwise information correlation: unconstrained ' ...
         '%.3f, at one variogram range %.3f, at one information horizon %.3f\n'], ...
        Des.info_corr_free_max, Des.info_corr_max, Des.info_corr_hor_max);

% ------------------------------------------ concentration versus pump test -
% A pumping test constrains the porosity at its own cell through the
% published log10 K to porosity relation used for the prior mean.  Its Fisher
% contribution is rank one in the KL basis with variance sigma_nK^2.
sd_lK = std(log10(P.Kw));
sd_nK = P.ne_network_sd * (cfg.sigma_K_rel / log(10)) / sd_lK;
Des.sd_pump_ne = sd_nK;
gain_pump = zeros(na, 1);
gain_conc = zeros(na, 1);
SigC = S.Sigma * c_goal;
for a = 1:na
    Ja = S.J(a, :);
    q = Ja * SigC;
    d = 1 + (Ja * S.Sigma * Ja') / sig2;
    gain_conc(a) = (q^2 / sig2) / d;
    pa = P.Phi(a, :);
    qp = pa * SigC;
    dp = 1 + (pa * S.Sigma * pa') / sd_nK^2;
    gain_pump(a) = (qp^2 / sd_nK^2) / dp;
end
Des.gain_conc = gain_conc;
Des.gain_pump = gain_pump;
Des.gain_conc_best = max(gain_conc(cand));
Des.gain_pump_best = max(gain_pump(cand));
Des.horizon_frac = cfg.horizon_frac;
Des.R_equiv = Des.gain_conc_best / Des.gain_pump_best;
fprintf(['[design] best single observation reduces Var(g) by %.4g d^2 ' ...
         '(concentration) vs %.4g d^2 (pumping test)\n'], ...
        max(gain_conc(cand)), max(gain_pump(cand)));
fprintf(['[design] information equivalence factor R = %.3f ' ...
         '(concentration observations per pumping test)\n'], Des.R_equiv);
end

% =========================================================================
function r = local_greedy(Sigma, J, c, cand, sig2, nsel, mode, min_sep, G)
sel = [];
logdet_gain = zeros(nsel, 1);
var_goal = zeros(nsel + 1, 1);
var_goal(1) = c' * Sigma * c;
Sg = Sigma;
for j = 1:nsel
    SgC = Sg * c;
    Jc = J * Sg;                      % na x m
    quad = sum(Jc .* J, 2);           % J Sg J'
    switch mode
        case 'D'
            score = log(1 + quad / sig2);
        case 'goal'
            num = (J * SgC) .^ 2 / sig2;
            score = num ./ (1 + quad / sig2);
    end
    score(~cand) = -inf;
    [~, best] = max(score);
    sel(end + 1, 1) = best; %#ok<AGROW>
    cand(best) = false;
    if min_sep > 0
        cand(hypot(G.xa - G.xa(best), G.ya - G.ya(best)) < min_sep) = false;
    end
    Jb = J(best, :);
    u = Sg * Jb';
    denom = sig2 + Jb * u;
    logdet_gain(j) = log(1 + (Jb * u) / sig2);
    Sg = Sg - (u * u') / denom;
    Sg = (Sg + Sg') / 2;
    var_goal(j + 1) = c' * Sg * c;
end
r.cells = sel;
r.logdet_gain = logdet_gain;
r.var_goal = var_goal;
r.Sigma_after = Sg;
end
