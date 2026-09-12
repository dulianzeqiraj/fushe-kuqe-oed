function M = mrva_replicate(cfg)
%MRVA_REPLICATE  Independent replication on a public USGS dataset.
%
%  Data: 44 wells of the Mississippi Alluvial Plain with hydraulic
%  conductivity and transmissivity from slug tests (Pugh, USGS Scientific
%  Investigations Report 2023-5101, doi:10.3133/sir20235101), together with
%  the DRASTIC ratings assigned to each well in the project file.  The
%  DRASTIC index of every well is recomputed here from its seven ratings and
%  the standard weights, and the recomputation must match the stored value
%  exactly or the run stops.
%
%  What is replicated is the design comparison, not the transport model.  No
%  flow field, boundary or nitrate record is available for this site, and
%  none is invented.  Both sites are therefore compared with the same local
%  observation operator, in which a new characterisation point measures the
%  field at its own location.  Every reported quantity is dimensionless or a
%  distance, so nothing depends on the absolute porosity scale, which this
%  dataset does not constrain.

T = readtable(fullfile(cfg.dir_data, 'mrva_44.csv'));
M.n_wells = height(T);

% ---------------------------------------------- recompute DRASTIC indices --
W = struct('D', 5, 'R', 4, 'A', 3, 'S', 2, 'T', 1, 'I', 5, 'C', 3);
V = W.D * T.D_rating + W.R * T.R_rating + W.A * T.A_rating + ...
    W.S * T.S_rating + W.T * T.T_rating + W.I * T.I_rating + W.C * T.C_rating;
M.drastic_mismatch = sum(abs(V - T.V_DRASTIC) > 1e-9);
if M.drastic_mismatch > 0
    error('mrva: %d DRASTIC indices do not reproduce from their ratings', ...
          M.drastic_mismatch);
end
fprintf('[mrva] %d wells, all DRASTIC indices reproduce from their ratings\n', ...
        M.n_wells);

x = T.x_km * 1000;  y = T.y_km * 1000;  K = T.K_m_per_d;
M.K_range = [min(K) max(K)];  M.K_mean = mean(K);  M.K_sd = std(K);
M.V_range = [min(V) max(V)];  M.V_mean = mean(V);  M.V_sd = std(V);
fprintf('[mrva] K %.2f - %.1f m/day (mean %.1f, sd %.1f); DRASTIC %g - %g (mean %.1f)\n', ...
        M.K_range(1), M.K_range(2), M.K_mean, M.K_sd, ...
        M.V_range(1), M.V_range(2), M.V_mean);

% ------------------------------------------------------------- variogram ---
lk = log10(K);
nb = 12;
[I, J] = find(triu(ones(numel(x)), 1));
h = hypot(x(I) - x(J), y(I) - y(J));
g = 0.5 * (lk(I) - lk(J)) .^ 2;
edges = linspace(0, 0.6 * max(h), nb + 1);
hb = zeros(nb, 1);  gb = nan(nb, 1);  nn = zeros(nb, 1);
for k = 1:nb
    s = h >= edges(k) & h < edges(k + 1);
    nn(k) = sum(s);
    if nn(k) > 0, hb(k) = mean(h(s));  gb(k) = mean(g(s)); end
end
ok = nn >= 10 & isfinite(gb);
s0 = var(lk);
obj = @(p) sum(nn(ok) .* (gb(ok) - (abs(p(1)) + abs(p(2)) * ...
              (1 - exp(-hb(ok) / abs(p(3)))))) .^ 2);
p = fminsearch(obj, [0.1 * s0, s0, 0.3 * max(hb)], ...
               optimset('Display', 'off', 'MaxFunEvals', 5000));
M.nugget = abs(p(1));  M.sill = abs(p(2));  M.range = abs(p(3));
fprintf('[mrva] exponential variogram on log10 K: nugget %.4f, sill %.4f, range %.0f m\n', ...
        M.nugget, M.sill, M.range);

% ------------------------------------------------------------------ grid ---
hull = convhull(x, y);
dx = cfg.mrva_dx;
xs = min(x):dx:max(x);
ys = min(y):dx:max(y);
[XX, YY] = meshgrid(xs, ys);
in = inpolygon(XX, YY, x(hull), y(hull));
M.xa = XX(in);  M.ya = YY(in);  M.n_active = nnz(in);
M.area_km2 = M.n_active * dx^2 / 1e6;
fprintf('[mrva] convex hull of the wells: %d cells of %g m, %.0f km2\n', ...
        M.n_active, dx, M.area_km2);

% -------------------------------------------------- vulnerability field ----
M.V = local_krige(x, y, V, M.xa, M.ya, M.nugget, M.sill, M.range);
fprintf('[mrva] kriged DRASTIC %.1f - %.1f (mean %.1f)\n', ...
        min(M.V), max(M.V), mean(M.V));

% ------------------------------------------------------------- KL basis ----
% Unit prior standard deviation: every reported quantity is a ratio or a
% distance, so the absolute scale cancels.
land = 1:cfg.nystrom_stride:M.n_active;
Xl = M.xa(land);  Yl = M.ya(land);
Cll = exp(-hypot(Xl - Xl', Yl - Yl') / M.range);
Cll = (Cll + Cll') / 2;
[U, Lam] = eig(Cll, 'vector');
[Lam, o] = sort(Lam, 'descend');  U = U(:, o);
m = min(cfg.n_kl, sum(Lam > 0));
Lam = Lam(1:m);  U = U(:, 1:m);
Phi = zeros(M.n_active, m);
for a = 1:20000:M.n_active
    b = min(a + 19999, M.n_active);
    Cxl = exp(-hypot(M.xa(a:b) - Xl', M.ya(a:b) - Yl') / M.range);
    Phi(a:b, :) = Cxl * (U ./ Lam');
end
for j = 1:m
    nr = norm(Phi(:, j));
    if nr > 0, Phi(:, j) = Phi(:, j) / nr; end
end
M.Phi = Phi;
M.lam = Lam * (M.n_active / numel(land));
M.kl_var_captured = sum(M.lam) / M.n_active;
fprintf('[mrva] KL basis: %d modes, %.1f%% of the prior variance\n', ...
        m, 100 * M.kl_var_captured);

% ---------------------------------------------------------------- design ---
cand = true(M.n_active, 1);
for k = 1:numel(x)
    [~, ci] = min(hypot(M.xa - x(k), M.ya - y(k)));
    cand(ci) = false;
end
w = (M.V / mean(M.V)) .^ cfg.gamma_Vdesign;
M.oed = oed_local(M.Phi, M.lam, w, cfg.sd_local, cand, cfg.n_new_wells, ...
                  M.xa, M.ya, M.range);
fprintf(['[mrva] variance of the vulnerability-weighted mean reduced by ' ...
         '%.1f%% (D-optimal) and %.1f%% (goal-oriented)\n'], ...
        100 * M.oed.var_reduction_D, 100 * M.oed.var_reduction_goal);
fprintf(['[mrva] goal-oriented design leaves %.1f%% less residual variance ' ...
         'than D-optimal; median separation %.0f m\n'], ...
        100 * M.oed.goal_advantage, M.oed.median_separation_m);
M.V_D = M.V(M.oed.D.cells);
M.V_goal = M.V(M.oed.goal.cells);
fprintf('[mrva] mean DRASTIC at chosen sites: D-optimal %.1f, goal-oriented %.1f\n', ...
        mean(M.V_D), mean(M.V_goal));
end

% =========================================================================
function zg = local_krige(x, y, z, xg, yg, nugget, sill, range)
n = numel(x);
Dd = hypot(x - x', y - y');
Gm = nugget + sill * (1 - exp(-Dd / range));
Gm(1:n+1:end) = 0;
A = [Gm, ones(n, 1); ones(1, n), 0] + 1e-10 * eye(n + 1);
zg = zeros(numel(xg), 1);
for a = 1:5000:numel(xg)
    b = min(a + 4999, numel(xg));
    Dg = hypot(xg(a:b) - x', yg(a:b) - y');
    Gg = nugget + sill * (1 - exp(-Dg / range));
    w = A \ [Gg'; ones(1, b - a + 1)];
    zg(a:b) = w(1:n, :)' * z;
end
end
