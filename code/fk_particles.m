function Par = fk_particles(cfg, G, P, F, Cal, S)
%FK_PARTICLES  Posterior characterisation by implicit sampling, against SIR.
%
%  The map from the porosity coefficients theta to the nitrate field is not
%  linear: theta enters through the storage term, so it changes how far the
%  front has travelled, and the observation operator is a nonlinear function
%  of it.  The posterior is therefore not Gaussian and a sampler is needed.
%
%  Cost function
%      F(theta) = 0.5 theta' diag(1/lambda) theta
%               + 0.5 || y - h(theta) ||^2 / sigma_obs^2
%
%  Implicit sampling (Chorin & Tu 2009; Chorin, Morzfeld & Tu 2010; Morzfeld
%  et al. 2012) minimises F, builds the quadratic F0 around the minimiser from
%  the Gauss-Newton Hessian, maps standard normal draws through the Cholesky
%  factor and corrects with the weight exp(-(F - F0)).  The comparison is with
%  sequential importance resampling drawing from the prior.
%
%  Two things are measured and reported: the effective sample size of each
%  sampler, and the number of iterations the minimisation takes under the
%  natural gradient of Amari (1998) versus the plain Euclidean gradient.

rng(cfg.seed);
m = numel(P.lam);
y = Cal.y_obs;
sig2 = cfg.sigma_obs ^ 2;
dt = cfg.cal_dt;
nstep = round(Cal.T0 / dt);

% forward evaluation at the observation wells for a given theta
    function c = hfun(th)
        ne = P.ne + P.Phi * th;
        ne = max(ne, 0.05);
        [A, Mass] = F.build(ne);
        Lh = decomposition(Mass / dt + A, 'lu');
        Md = Mass / dt;
        C = zeros(G.n_active, 1);
        for it = 1:nstep
            C = Lh \ (Md * C + Cal.source);
        end
        c = C(Cal.cell_of);
    end
    function v = Ffun(th)
        r = y - hfun(th);
        v = 0.5 * sum((th .^ 2) ./ P.lam) + 0.5 * sum(r .^ 2) / sig2;
    end

nfwd = 0;

% ------------------------------------------------------ Gauss-Newton MAP ---
% The Jacobian is held at its prior-state value, which is the Gauss-Newton
% approximation used throughout; the residual is recomputed each step.
th = zeros(m, 1);
Hgn = S.Itot;
Rch = chol(Hgn, 'lower');
Par.F_path = zeros(0, 1);
for k = 1:cfg.gn_iters
    r = y - hfun(th);  nfwd = nfwd + 1;
    g = th ./ P.lam - (S.Jobs' * r) / sig2;
    dth = -(Rch' \ (Rch \ g));
    % backtracking on the true cost
    a = 1;  f0 = Ffun(th);  nfwd = nfwd + 1;
    for bt = 1:8
        if Ffun(th + a * dth) < f0, break; end
        nfwd = nfwd + 1;
        a = a / 2;
    end
    th = th + a * dth;
    Par.F_path(end + 1, 1) = Ffun(th);  nfwd = nfwd + 1;
    if norm(dth) < 1e-10 * norm(th) + 1e-12, break; end
end
Par.theta_map = th;
Par.phi_F = Ffun(th);  nfwd = nfwd + 1;
Par.gn_iters_used = numel(Par.F_path);
fprintf('[part] Gauss-Newton MAP reached in %d iterations, F = %.4f\n', ...
        Par.gn_iters_used, Par.phi_F);

% --------------------------------- natural gradient versus Euclidean -------
% Both descend the same quadratic model of F around the prior state, so the
% comparison isolates the effect of the metric.  The Euclidean step uses the
% largest stable constant step, 1/L with L the largest eigenvalue.
Hs = (S.Itot + S.Itot') / 2;
ev = eig(Hs);
Par.cond_H = max(ev) / min(ev);
g0 = -(S.Jobs' * (y - hfun(zeros(m, 1)))) / sig2;  nfwd = nfwd + 1;
tol = 1e-6;
% natural gradient: preconditioned by the inverse Fisher information
xg = zeros(m, 1);  it_nat = 0;
while it_nat < cfg.max_desc
    gk = Hs * xg + g0;
    if sqrt(gk' * (Hs \ gk)) < tol * sqrt(g0' * (Hs \ g0)), break; end
    xg = xg - (Hs \ gk);
    it_nat = it_nat + 1;
end
% Euclidean gradient
xe = zeros(m, 1);  it_euc = 0;  step = 1 / max(ev);
while it_euc < cfg.max_desc
    gk = Hs * xe + g0;
    if sqrt(gk' * (Hs \ gk)) < tol * sqrt(g0' * (Hs \ g0)), break; end
    xe = xe - step * gk;
    it_euc = it_euc + 1;
end
Par.iters_natural = it_nat;
Par.iters_euclidean = it_euc;
Par.speedup = it_euc / max(it_nat, 1);
fprintf(['[part] condition number of the Gauss-Newton Hessian %.3e\n'], Par.cond_H);
fprintf(['[part] descent to a relative Riemannian gradient of %.0e: ' ...
         'natural %d iterations, Euclidean %d, ratio %.1f (sqrt(kappa) = %.1f)\n'], ...
        tol, it_nat, it_euc, Par.speedup, sqrt(Par.cond_H));

% ------------------------------------------------------ implicit sampling --
L = chol(Hgn, 'lower');
M = cfg.M_particles;
Xi = randn(m, M);
W = zeros(M, 1);
Fv = zeros(M, 1);  F0v = zeros(M, 1);
for k = 1:M
    xk = Par.theta_map + L' \ Xi(:, k);
    Fv(k) = Ffun(xk);  nfwd = nfwd + 1;
    dx = xk - Par.theta_map;
    F0v(k) = Par.phi_F + 0.5 * (dx' * (Hgn * dx));
    W(k) = exp(-(Fv(k) - F0v(k)));
end
Par.w_implicit = W;
Par.ess_implicit = (sum(W) ^ 2) / (M * sum(W .^ 2));
fprintf('[part] implicit sampling with M = %d: normalised ESS %.3f\n', ...
        M, Par.ess_implicit);

% posterior mean and spread of the goal-relevant field
Xs = Par.theta_map + (L' \ Xi);
wn = W / sum(W);
Par.theta_mean = Xs * wn;
dev = Xs - Par.theta_mean;
Par.theta_cov_diag = (dev .^ 2) * wn;

% --------------------------------------------------------------- SIR -------
Ms = cfg.M_sir;
Ws = zeros(Ms, 1);
for k = 1:Ms
    xk = sqrt(P.lam) .* randn(m, 1);
    r = y - hfun(xk);  nfwd = nfwd + 1;
    Ws(k) = exp(-0.5 * sum(r .^ 2) / sig2);
end
Ws = Ws / max(Ws);
Par.ess_sir = (sum(Ws) ^ 2) / (Ms * sum(Ws .^ 2));
fprintf('[part] SIR from the prior with M = %d: normalised ESS %.3f\n', ...
        Ms, Par.ess_sir);
Par.n_forward = nfwd;
fprintf('[part] %d forward transport solves in total\n', nfwd);
end
