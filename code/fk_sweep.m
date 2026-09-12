function Sw = fk_sweep(cfg, D, G, P, F, Cal, S, Des)
%FK_SWEEP  Do the conclusions survive the assumptions they rest on?
%
%  Three quantities are followed: the fraction of the total Fisher information
%  contributed by the 24 nitrate observations, the information-equivalence
%  factor between a concentration observation and a pumping test, and the
%  separation between the classical and goal-oriented designs.  Each is
%  recomputed as the assumptions marked [ASSUM] in fk_config are varied one at
%  a time.  The dispersivity and the flow azimuth require a new Jacobian and
%  are therefore the expensive ones.

Sw = struct();
sig2 = cfg.sigma_obs ^ 2;

% ------------------------------------ measurement error of the nitrate data
Sw.sigma_obs = [0.5 1.0 1.5 2.0];
Sw.frac_obs_vs_sigma = zeros(size(Sw.sigma_obs));
Sw.R_vs_sigma = zeros(size(Sw.sigma_obs));
for k = 1:numel(Sw.sigma_obs)
    s2 = Sw.sigma_obs(k) ^ 2;
    Iobs = (S.Jobs' * S.Jobs) / s2;
    Sw.frac_obs_vs_sigma(k) = trace(Iobs) / (trace(S.Iprior) + trace(Iobs));
    Sg = inv(S.Iprior + Iobs);
    Sg = (Sg + Sg') / 2;
    Sw.R_vs_sigma(k) = local_R(cfg, P, S, Des, Sg, s2);
    fprintf(['[sweep] sigma_obs %.1f mg/L : observation share of the Fisher ' ...
             'trace %.2e, R = %.3f\n'], Sw.sigma_obs(k), ...
            Sw.frac_obs_vs_sigma(k), Sw.R_vs_sigma(k));
end

% ---------------------------------------- relative error of a pumping test
Sw.sigma_K_rel = [0.05 0.10 0.15 0.20];
Sw.R_vs_sigmaK = zeros(size(Sw.sigma_K_rel));
for k = 1:numel(Sw.sigma_K_rel)
    c2 = cfg;  c2.sigma_K_rel = Sw.sigma_K_rel(k);
    Sw.R_vs_sigmaK(k) = local_R(c2, P, S, Des, S.Sigma, sig2);
    fprintf('[sweep] pumping-test relative error %.2f : R = %.3f\n', ...
            Sw.sigma_K_rel(k), Sw.R_vs_sigmaK(k));
end

% ------------------------------- weighting exponent of the goal functional
Sw.gammaV = [0 0.5 1.0 1.5 2.0 3.0];
Sw.sep_vs_gamma = zeros(size(Sw.gammaV));
Sw.adv_vs_gamma = zeros(size(Sw.gammaV));
Sw.V_goal_vs_gamma = zeros(size(Sw.gammaV));
xi_a = G.xa * G.exi(1) + G.ya * G.exi(2);
Lrem = max(xi_a) - xi_a;
cand = true(G.n_active, 1);  cand(Cal.cell_of) = false;
for k = 1:numel(Sw.gammaV)
    w = (P.V / mean(P.V)) .^ Sw.gammaV(k);
    c = (P.Phi' * (w .* Lrem)) / (F.q0 * sum(w));
    rD = local_greedy(S.Sigma, S.J, c, cand, sig2, cfg.n_new_wells, 'D', P.vrange, G);
    rG = local_greedy(S.Sigma, S.J, c, cand, sig2, cfg.n_new_wells, 'goal', P.vrange, G);
    d = zeros(cfg.n_new_wells, 1);
    for j = 1:cfg.n_new_wells
        d(j) = min(hypot(G.xa(rG.cells(j)) - G.xa(rD.cells), ...
                         G.ya(rG.cells(j)) - G.ya(rD.cells)));
    end
    Sw.sep_vs_gamma(k) = median(d);
    Sw.adv_vs_gamma(k) = (rD.var_goal(end) - rG.var_goal(end)) / rD.var_goal(end);
    Sw.V_goal_vs_gamma(k) = mean(P.V(rG.cells));
    fprintf(['[sweep] gamma_V %.1f : median separation %.0f m, goal-oriented ' ...
             'advantage %.1f%%, mean DRASTIC at sites %.1f\n'], ...
            Sw.gammaV(k), Sw.sep_vs_gamma(k), 100 * Sw.adv_vs_gamma(k), ...
            Sw.V_goal_vs_gamma(k));
end

% ------------------------------------------------------------ dispersivity
Sw.alphaL = [100 250 500 1000];
Sw.frac_obs_vs_alphaL = zeros(size(Sw.alphaL));
Sw.R_vs_alphaL = zeros(size(Sw.alphaL));
Sw.sep_vs_alphaL = zeros(size(Sw.alphaL));
for k = 1:numel(Sw.alphaL)
    c2 = cfg;  c2.alphaL = Sw.alphaL(k);
    F2 = fk_forward_quiet(c2, G, P);
    S2 = fk_sensitivity_quiet(c2, G, P, F2, Cal);
    Sw.frac_obs_vs_alphaL(k) = S2.frac_obs;
    Sw.R_vs_alphaL(k) = local_R(c2, P, S2, Des, S2.Sigma, c2.sigma_obs^2);
    w = (P.V / mean(P.V)) .^ cfg.gamma_Vdesign;
    c = (P.Phi' * (w .* Lrem)) / (F2.q0 * sum(w));
    rD = local_greedy(S2.Sigma, S2.J, c, cand, c2.sigma_obs^2, cfg.n_new_wells, 'D', P.vrange, G);
    rG = local_greedy(S2.Sigma, S2.J, c, cand, c2.sigma_obs^2, cfg.n_new_wells, 'goal', P.vrange, G);
    d = zeros(cfg.n_new_wells, 1);
    for j = 1:cfg.n_new_wells
        d(j) = min(hypot(G.xa(rG.cells(j)) - G.xa(rD.cells), ...
                         G.ya(rG.cells(j)) - G.ya(rD.cells)));
    end
    Sw.sep_vs_alphaL(k) = median(d);
    fprintf(['[sweep] alphaL %4d m : observation share %.2e, R = %.3f, ' ...
             'median separation %.0f m\n'], Sw.alphaL(k), ...
            Sw.frac_obs_vs_alphaL(k), Sw.R_vs_alphaL(k), Sw.sep_vs_alphaL(k));
end
end

% =========================================================================
function R = local_R(cfg, P, S, Des, Sigma, sig2)
sd_lK = std(log10(P.Kw));
sd_nK = P.ne_network_sd * (cfg.sigma_K_rel / log(10)) / sd_lK;
c = Des.c_goal;
SigC = Sigma * c;
q  = S.J * SigC;
qd = sum((S.J * Sigma) .* S.J, 2);
gc = (q .^ 2 / sig2) ./ (1 + qd / sig2);
qp  = P.Phi * SigC;
qpd = sum((P.Phi * Sigma) .* P.Phi, 2);
gp = (qp .^ 2 / sd_nK^2) ./ (1 + qpd / sd_nK^2);
R = max(gc) / max(gp);
end

% =========================================================================
function r = local_greedy(Sigma, J, c, cand, sig2, nsel, mode, min_sep, G)
sel = zeros(nsel, 1);
var_goal = zeros(nsel + 1, 1);
var_goal(1) = c' * Sigma * c;
Sg = Sigma;
for j = 1:nsel
    SgC = Sg * c;
    quad = sum((J * Sg) .* J, 2);
    switch mode
        case 'D',    score = log(1 + quad / sig2);
        case 'goal', score = ((J * SgC) .^ 2 / sig2) ./ (1 + quad / sig2);
    end
    score(~cand) = -inf;
    [~, b] = max(score);
    sel(j) = b;  cand(b) = false;
    if min_sep > 0
        cand(hypot(G.xa - G.xa(b), G.ya - G.ya(b)) < min_sep) = false;
    end
    Jb = J(b, :);  u = Sg * Jb';
    Sg = Sg - (u * u') / (sig2 + Jb * u);
    Sg = (Sg + Sg') / 2;
    var_goal(j + 1) = c' * Sg * c;
end
r.cells = sel;  r.var_goal = var_goal;
end

% =========================================================================
function F = fk_forward_quiet(cfg, G, P)
ev = evalc('F = fk_forward(cfg, G, P);'); %#ok<NASGU>
end

function S = fk_sensitivity_quiet(cfg, G, P, F, Cal)
ev = evalc('S = fk_sensitivity(cfg, G, P, F, Cal);'); %#ok<NASGU>
end
