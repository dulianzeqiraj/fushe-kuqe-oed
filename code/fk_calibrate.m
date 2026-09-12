function Cal = fk_calibrate(cfg, D, G, P, F)
%FK_CALIBRATE  Estimate the nitrate loading field and its onset time.
%
%  A first attempt used a single loading amplitude with the spatial shape of
%  the DRASTIC index.  It failed: the best fit over the eight exponents in
%  cfg.cal_kappas reached R2 = 0.04 and a negative leave-one-out R2, and the
%  onset time collapsed onto the shortest step in the scan.  The reason is
%  visible in the data.  Nitrate in this aquifer is organised north to south
%  (correlation with northing -0.57 across the 31 instrumented wells) far more
%  strongly than it is organised by vulnerability (correlation with the
%  DRASTIC index +0.30).  The pattern is set by where the nitrogen enters, not
%  by how the aquifer transports it, so a loading field whose shape is fixed a
%  priori cannot reproduce it.
%
%  The loading is therefore expanded in a small set of large-scale spatial
%  functions and its coefficients estimated from the observations:
%
%     m_s(x) = sum_j beta_j psi_j(x),
%     psi = { 1, xi, eta, xi^2, eta^2, xi*eta, V }
%
%  with xi and eta the along- and across-flow coordinates normalised to
%  [-1, 1] and V the normalised DRASTIC index.  The transport equation is
%  linear in the source, so the response to each psi_j is computed once and
%  the coefficients follow from linear least squares at every candidate onset
%  time.  The DRASTIC coefficient is then a result rather than an assumption:
%  it answers whether vulnerability explains any of the loading once the
%  regional trend is removed.

rng(cfg.seed);

% ------------------------------------------------- observation locations ---
ob = D.inside & D.has_NO3;
Cal.obs_id = D.well_id(ob);
xo = D.x(ob);  yo = D.y(ob);  Cal.y_obs = D.NO3(ob);
Cal.n_obs = numel(Cal.y_obs);
cell_of = zeros(Cal.n_obs, 1);
for k = 1:Cal.n_obs
    [~, cell_of(k)] = min(hypot(G.xa - xo(k), G.ya - yo(k)));
end
Cal.cell_of = cell_of;
Cal.x_obs = xo;  Cal.y_obs_coord = yo;
fprintf('[cal] %d observation wells inside the domain, max well-to-cell %.0f m\n', ...
        Cal.n_obs, max(hypot(G.xa(cell_of) - xo, G.ya(cell_of) - yo)));

% ------------------------------- the null model that failed, run and reported
% A loading whose spatial shape is fixed to a power of the DRASTIC index, with
% one amplitude and one onset time.  It is kept in the pipeline because the
% paper reports that it fails; the numbers below are that failure.
[A0, M0] = F.build(P.ne);
Lhs0 = decomposition(M0 / cfg.cal_dt + A0, 'lu');
Mdt0 = M0 / cfg.cal_dt;
tg0 = (1:cfg.cal_nsteps)' * cfg.cal_dt;
Vr = P.V / mean(P.V);
sst0 = sum((Cal.y_obs - mean(Cal.y_obs)) .^ 2);
Cal.null_kappa = cfg.cal_kappas;
Cal.null_R2 = zeros(size(cfg.cal_kappas));
Cal.null_T0 = zeros(size(cfg.cal_kappas));
for kk = 1:numel(cfg.cal_kappas)
    s0 = (cfg.dx * cfg.dy) * (Vr .^ cfg.cal_kappas(kk));
    C0 = zeros(G.n_active, 1);
    H0 = zeros(cfg.cal_nsteps, Cal.n_obs);
    for it = 1:cfg.cal_nsteps
        C0 = Lhs0 \ (Mdt0 * C0 + s0);
        H0(it, :) = C0(cell_of)';
    end
    num0 = H0 * Cal.y_obs;  den0 = sum(H0 .^ 2, 2);
    sse0 = sum(Cal.y_obs .^ 2) - (num0 .^ 2) ./ max(den0, realmin);
    [b0, i0] = min(sse0);
    Cal.null_R2(kk) = 1 - b0 / sst0;
    Cal.null_T0(kk) = tg0(i0);
end
[Cal.null_best_R2, ib] = max(Cal.null_R2);
Cal.null_best_kappa = cfg.cal_kappas(ib);
Cal.null_best_T0 = Cal.null_T0(ib);
fprintf(['[cal] null model with the loading shaped as DRASTIC^kappa: best ' ...
         'R2 %.3f at kappa %.2f, T0 %.0f d\n'], Cal.null_best_R2, ...
        Cal.null_best_kappa, Cal.null_best_T0);

% ------------------------------------------------------ source basis -------
xi_a  = G.xa * G.exi(1) + G.ya * G.exi(2);
eta_a = G.xa * G.eta(1) + G.ya * G.eta(2);
nz = @(v) 2 * (v - min(v)) / (max(v) - min(v)) - 1;
xs = nz(xi_a);  es = nz(eta_a);  vs = nz(P.V);
Psi = [ones(G.n_active, 1), xs, es, xs.^2, es.^2, xs.*es, vs];
Cal.psi_names = {'1', 'xi', 'eta', 'xi^2', 'eta^2', 'xi*eta', 'V_DRASTIC'};
nb = size(Psi, 2);

% ------------------------------------------------------------ marching -----
[A, Mass] = F.build(P.ne);
dt = cfg.cal_dt;  nstep = cfg.cal_nsteps;
Lhs = decomposition(Mass / dt + A, 'lu');
Mdt = Mass / dt;
tgrid = (1:nstep)' * dt;

Hobs = zeros(nstep, Cal.n_obs, nb);       % response at the wells per unit beta
tS = tic;
for j = 1:nb
    s = (cfg.dx * cfg.dy) * Psi(:, j);
    C = zeros(G.n_active, 1);
    for it = 1:nstep
        C = Lhs \ (Mdt * C + s);
        Hobs(it, :, j) = C(cell_of)';
    end
end
fprintf('[cal] %d source-basis marches of %d steps in %.0f s\n', nb, nstep, toc(tS));

% ---------------------------------------- least squares at every onset time
y = Cal.y_obs;
sst = sum((y - mean(y)) .^ 2);
sse_t = inf(nstep, 1);
beta_t = zeros(nb, nstep);
for it = 1:nstep
    X = squeeze(Hobs(it, :, :));
    if rcond(X' * X) < 1e-14, continue; end
    b = (X' * X) \ (X' * y);
    beta_t(:, it) = b;
    sse_t(it) = sum((y - X * b) .^ 2);
end
[sse, it_best] = min(sse_t);
Cal.T0 = tgrid(it_best);
Cal.beta = beta_t(:, it_best);
Cal.rmse = sqrt(sse / Cal.n_obs);
Cal.R2 = 1 - sse / sst;
Cal.n_par = nb + 1;                        % the basis plus the onset time
Cal.R2_adj = 1 - (1 - Cal.R2) * (Cal.n_obs - 1) / (Cal.n_obs - Cal.n_par - 1);
Cal.sse_profile = sse_t;
Cal.tgrid = tgrid;

fprintf(['[cal] onset time T0 = %.0f d (%.1f yr), RMSE %.3f mg/L, ' ...
         'R2 %.3f, adjusted R2 %.3f\n'], Cal.T0, Cal.T0 / 365.25, ...
        Cal.rmse, Cal.R2, Cal.R2_adj);
for j = 1:nb
    fprintf('[cal]   beta(%-9s) = %+.4e\n', Cal.psi_names{j}, Cal.beta(j));
end

% ------------------------------------------------------- source field ------
Cal.source_density = Psi * Cal.beta;                    % mg/L/day
Cal.frac_negative = mean(Cal.source_density < 0);
fprintf('[cal] fitted loading %.3e to %.3e mg/L/d, negative over %.1f%% of the domain\n', ...
        min(Cal.source_density), max(Cal.source_density), ...
        100 * Cal.frac_negative);
% A loading field must not be a sink.  Where least squares drives it below
% zero it is clipped and the amplitude refitted on the clipped shape.
shape = max(Cal.source_density, 0);
s = (cfg.dx * cfg.dy) * shape;
C = zeros(G.n_active, 1);
Hc = zeros(nstep, Cal.n_obs);
for it = 1:nstep
    C = Lhs \ (Mdt * C + s);
    Hc(it, :) = C(cell_of)';
end
num = Hc * y;  den = sum(Hc .^ 2, 2);
a_t = num ./ max(den, realmin);
sse_c = sum(y .^ 2) - (num .^ 2) ./ max(den, realmin);
[sse2, it2] = min(sse_c);
Cal.T0_clipped = tgrid(it2);
Cal.amp_clipped = a_t(it2);
Cal.rmse_clipped = sqrt(sse2 / Cal.n_obs);
Cal.R2_clipped = 1 - sse2 / sst;
fprintf(['[cal] non-negative loading: T0 %.0f d (%.1f yr), RMSE %.3f mg/L, ' ...
         'R2 %.3f\n'], Cal.T0_clipped, Cal.T0_clipped / 365.25, ...
        Cal.rmse_clipped, Cal.R2_clipped);

Cal.source = (cfg.dx * cfg.dy) * Cal.amp_clipped * shape;
Cal.T0 = Cal.T0_clipped;
Cal.rmse = Cal.rmse_clipped;
Cal.R2 = Cal.R2_clipped;

C = zeros(G.n_active, 1);
for it = 1:it2
    C = Lhs \ (Mdt * C + Cal.source);
end
Cal.C = C;
Cal.y_fit = C(cell_of);
Cal.resid = y - Cal.y_fit;
fprintf('[cal] modelled nitrate %.2f - %.2f mg/L (observed %.2f - %.2f)\n', ...
        min(C), max(C), min(y), max(y));

% ---------------------------------------------- what DRASTIC alone explains
% Same fit with the vulnerability column removed, to isolate its contribution.
keep = 1:(nb - 1);
sse_t2 = inf(nstep, 1);
for it = 1:nstep
    X = squeeze(Hobs(it, :, keep));
    if rcond(X' * X) < 1e-14, continue; end
    b = (X' * X) \ (X' * y);
    sse_t2(it) = sum((y - X * b) .^ 2);
end
Cal.R2_noV = 1 - min(sse_t2) / sst;
fprintf(['[cal] dropping the DRASTIC column changes R2 from %.3f to %.3f ' ...
         '(partial contribution %.3f)\n'], ...
        1 - sse / sst, Cal.R2_noV, (1 - sse / sst) - Cal.R2_noV);

% ------------------------------------------------------------- LOO-CV ------
pred = zeros(Cal.n_obs, 1);
for k = 1:Cal.n_obs
    kp = true(Cal.n_obs, 1);  kp(k) = false;
    best = inf;  bp = 0;
    for it = 1:nstep
        X = squeeze(Hobs(it, kp, :));
        if rcond(X' * X) < 1e-14, continue; end
        b = (X' * X) \ (X' * y(kp));
        e = sum((y(kp) - X * b) .^ 2);
        if e < best
            best = e;
            bp = squeeze(Hobs(it, k, :))' * b;
        end
    end
    pred(k) = bp;
end
Cal.loo_pred = pred;
Cal.loo_rmse = sqrt(mean((y - pred) .^ 2));
Cal.loo_R2 = 1 - sum((y - pred) .^ 2) / sst;
fprintf('[cal] leave-one-out: RMSE %.3f mg/L, R2 %.3f\n', ...
        Cal.loo_rmse, Cal.loo_R2);
end
