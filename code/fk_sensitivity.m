function S = fk_sensitivity(cfg, G, P, F, Cal)
%FK_SENSITIVITY  Jacobian of the nitrate field with respect to the porosity
%                field, expressed in the Karhunen-Loeve basis of the prior.
%
%  The porosity field is written  n_e(x) = n_prior(x) + sum_k theta_k phi_k(x)
%  with theta ~ N(0, diag(lambda)).  Differentiating the discrete transport
%  equation with respect to theta_k gives
%
%    (Mass/dt + A) Ck^n = (Mass/dt) Ck^(n-1)
%                       - diag(phi_k * cellVolume) (C^n - C^(n-1)) / dt
%
%  which is the same left-hand side as the base march.  One sparse
%  factorisation therefore serves the base solve and all KL modes, and the
%  derivative is exact for the discrete system rather than a finite
%  difference.  A finite-difference check on a random mode is run and
%  reported.

na = G.n_active;
m = numel(P.lam);
dt = cfg.cal_dt;
nstep = round(Cal.T0 / dt);
fprintf('[sens] %d KL modes, %d time steps to T0 = %.0f d\n', m, nstep, Cal.T0);

[A, Mass] = F.build(P.ne);
Lhs = decomposition(Mass / dt + A, 'lu');
Mdt = Mass / dt;
s = Cal.source;

% ------------------------------------------------------------ base march ---
Cbase = zeros(na, nstep + 1);
C = zeros(na, 1);
for it = 1:nstep
    C = Lhs \ (Mdt * C + s);
    Cbase(:, it + 1) = C;
end
S.C = C;
fprintf('[sens] base field at T0: %.3f - %.3f mg/L\n', min(C), max(C));

% ------------------------------------------------------- mode sensitivity --
vol = F.cell_volume;
J = zeros(na, m);
tS = tic;
for k = 1:m
    wk = P.Phi(:, k) * vol / dt;         % diag(phi_k * V) / dt
    Ck = zeros(na, 1);
    for it = 1:nstep
        rhs = Mdt * Ck - wk .* (Cbase(:, it + 1) - Cbase(:, it));
        Ck = Lhs \ rhs;
    end
    J(:, k) = Ck;
    if mod(k, 25) == 0
        fprintf('[sens] mode %d/%d, %.0f s elapsed\n', k, m, toc(tS));
    end
end
S.J = J;
fprintf('[sens] Jacobian built in %.0f s\n', toc(tS));

% -------------------------------------------------- finite-difference check
% Run at three step sizes.  The derivative is exact for the discrete system, so
% the residual should fall linearly with the step, which is what distinguishes
% a correct derivative from a coincidence at one step.
kchk = min(3, m);
S.fd_steps = [0.25 0.05 0.01] * sqrt(P.lam(kchk));
S.fd_rel = zeros(size(S.fd_steps));
ad = J(:, kchk);
for q = 1:numel(S.fd_steps)
    eps_k = S.fd_steps(q);
    [A2, Mass2] = F.build(P.ne + eps_k * P.Phi(:, kchk));
    Lhs2 = decomposition(Mass2 / dt + A2, 'lu');
    Mdt2 = Mass2 / dt;
    C2 = zeros(na, 1);
    for it = 1:nstep
        C2 = Lhs2 \ (Mdt2 * C2 + s);
    end
    fd = (C2 - S.C) / eps_k;
    S.fd_rel(q) = norm(fd - ad) / max(norm(ad), realmin);
    fprintf(['[sens] finite-difference check on mode %d, step %.3e: ' ...
             'relative difference %.3e\n'], kchk, eps_k, S.fd_rel(q));
end
S.fd_check_rel = S.fd_rel(end);
S.fd_slope = log(S.fd_rel(1) / S.fd_rel(end)) / log(S.fd_steps(1) / S.fd_steps(end));
fprintf('[sens] residual falls with step as a power of %.2f (first order is 1)\n', ...
        S.fd_slope);

% ------------------------------------------------- Fisher information ------
S.Iprior = diag(1 ./ P.lam);
Jo = J(Cal.cell_of, :);
S.Jobs = Jo;
S.Iobs = (Jo' * Jo) / cfg.sigma_obs^2;
S.Itot = S.Iprior + S.Iobs;
S.Sigma = inv(S.Itot);
S.Sigma = (S.Sigma + S.Sigma') / 2;

S.tr_prior = trace(S.Iprior);
S.tr_obs   = trace(S.Iobs);
S.frac_prior = S.tr_prior / (S.tr_prior + S.tr_obs);
S.frac_obs   = 1 - S.frac_prior;

% Degrees of freedom for signal: how many directions of the parameter space
% the observations actually constrain.  This is the measure that should be
% quoted.  The trace ratio above is reported as well because earlier work on
% this aquifer uses it, but it is biased: the trace of the prior information
% is dominated by the smallest prior variances, that is by the directions the
% prior already pins down, so it flatters the prior by construction.  The
% degrees of freedom are scale free and bounded by the number of modes.
S.dofs = trace(S.Iobs * S.Sigma);
S.n_modes = m;
fprintf(['[sens] degrees of freedom for signal %.4f out of %d modes, ' ...
         'so the observations constrain the equivalent of %.4f directions\n'], ...
        S.dofs, m, S.dofs);
S.logdet_prior = sum(log(1 ./ P.lam));
S.logdet_tot   = 2 * sum(log(diag(chol(S.Itot))));
S.info_gain_nats = 0.5 * (S.logdet_tot - S.logdet_prior);
fprintf(['[sens] Fisher trace: prior %.3e (%.1f%%), observations %.3e ' ...
         '(%.1f%%)\n'], S.tr_prior, 100*S.frac_prior, S.tr_obs, 100*S.frac_obs);
fprintf(['[sens] log|I_total| - log|I_prior| = %.3f, so the 24 nitrate wells ' ...
         'carry %.3f nats\n'], S.logdet_tot - S.logdet_prior, S.info_gain_nats);

% Posterior sd of the porosity field, back in space.  The posterior covariance
% is not diagonal in the KL basis, so the pointwise variance is the diagonal of
% Phi * Sigma * Phi', computed row by row rather than from diag(Sigma) alone.
S.var_post = sum((P.Phi * S.Sigma) .* P.Phi, 2);
S.sd_post = sqrt(max(S.var_post, 0));
S.var_prior_field = sum((P.Phi .^ 2) .* (P.lam'), 2);
S.sd_prior_field = sqrt(max(S.var_prior_field, 0));
red = 1 - mean(S.sd_post) / mean(S.sd_prior_field);
fprintf(['[sens] mean prior sd in the KL subspace %.5f -> posterior %.5f ' ...
         '(%.1f%% reduction)\n'], mean(S.sd_prior_field), mean(S.sd_post), ...
        100 * red);
S.sd_reduction = red;
end
