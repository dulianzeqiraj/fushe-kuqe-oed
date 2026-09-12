function P = fk_prior_field(cfg, D, G)
%FK_PRIOR_FIELD  Porosity prior, vulnerability field and the KL basis.
%
%  The porosity prior is NOT re-derived here.  The Dual-Pathway fusion for this
%  aquifer is already published (Zeqiraj et al. 2026, J Hazard Mater Adv 23,
%  101261; Zeqiraj & Beqiraj 2026, J Contam Hydrol 283, 105086).  This routine
%  consumes that result: the facies reference ranges and the pathway-error
%  correlations are read from cfg, the measured K of each well places it inside
%  its facies range, and the correlation-aware variance of the cited work sets
%  the prior spread.
%
%  Outputs
%    P.ne_well   porosity prior at each of the 43 wells with a measured K
%    P.sd_well   prior sd at each of those wells (facies dependent)
%    P.sd_iv     the independence-assumption sd, for the comparison in the text
%    P.vrange    correlation length fitted to the measured log10 K, m
%    P.ne        kriged porosity prior on the grid
%    P.sd        kriged prior sd on the grid
%    P.V         kriged DRASTIC vulnerability on the grid
%    P.Phi,P.lam KL basis of the prior covariance (n_active x n_kl, n_kl x 1)

kk = D.has_K & D.facies > 0;
P.well_idx = find(kk);
nK = numel(P.well_idx);
xw = D.x(kk);  yw = D.y(kk);  Kw = D.K(kk);  fw = D.facies(kk);

% ---------------------------------------- does the lithology label carry K? --
% The porosity prior would normally be keyed on facies.  In this workbook the
% lithological descriptor and the pumping-test K are almost independent, so a
% facies-keyed prior would impose structure the data do not contain.  The test
% is reported because it governs the choice made immediately below.
Kf_mean = zeros(1, 3);  Kf_sd = zeros(1, 3);  Kf_n = zeros(1, 3);
for f = 1:3
    s = (fw == f);
    Kf_n(f) = sum(s);  Kf_mean(f) = mean(Kw(s));  Kf_sd(f) = std(Kw(s));
    fprintf('[prior] facies %-7s n=%2d  K %3g-%3g  mean %.1f sd %.1f m/day\n', ...
            cfg.facies_short{f}, Kf_n(f), min(Kw(s)), max(Kw(s)), ...
            Kf_mean(f), Kf_sd(f));
end
grand = mean(Kw);
ss_between = sum(Kf_n .* (Kf_mean - grand).^2);
ss_total   = sum((Kw - grand).^2);
P.eta2_facies_K = ss_between / ss_total;
dfb = 2;  dfw = nK - 3;
F = (ss_between / dfb) / ((ss_total - ss_between) / dfw);
P.F_facies_K = F;
P.p_facies_K = 1 - fcdf(F, dfb, dfw);
fprintf(['[prior] one-way ANOVA of K on the lithology label: ' ...
         'eta^2 = %.4f, F(%d,%d) = %.3f, p = %.3f\n'], ...
        P.eta2_facies_K, dfb, dfw, F, P.p_facies_K);

% ------------------------------------------------- prior mean at the wells --
% The label carries almost no information about K, so the spatial pattern of
% the prior is taken from the measurement that does: a monotone map of log10 K
% onto the published network porosity statistics.  Both constants come from
% the fusion already published for this aquifer (network mean 0.282, network
% range 0.22-0.35; Zeqiraj et al. 2026, J Hazard Mater Adv 23, 101261,
% Table 3).  Nothing here is fitted to any result of the present paper.
lK = log10(Kw);
z  = (lK - mean(lK)) / std(lK);
P.ne_network_mean = 0.282;            % [PUB]
P.ne_network_sd   = 0.020;            % [DERIV] published range 0.22-0.35 over
                                      %         a 43-well network is +-3.2 sd
ne_w = P.ne_network_mean + P.ne_network_sd * z;
ne_w = min(max(ne_w, 0.22), 0.35);    % [PUB] published physical range
fprintf('[prior] prior mean from log10 K: %.4f - %.4f, network mean %.4f\n', ...
        min(ne_w), max(ne_w), mean(ne_w));

% ------------------------------------------------------- prior sd ----------
% Independence (inverse-variance) combination, and the correlation-aware BLUE
% variance of Zeqiraj & Beqiraj (2026), Eq. (2).  Because the facies labels do
% not separate this dataset, the network-representative medium-facies
% correlation is used, and rho is swept over the full published range in the
% sensitivity analysis.
s1 = cfg.sigma_VS;  s2 = cfg.sigma_KC;
P.sd_iv = sqrt(s1^2 * s2^2 / (s1^2 + s2^2));
sd_facies = zeros(1, 3);
for f = 1:3
    r = cfg.rho_pathway(f);
    sd_facies(f) = sqrt((1 - r^2) * s1^2 * s2^2 / (s1^2 + s2^2 - 2*r*s1*s2));
end
P.rho_used = cfg.rho_pathway(2);
sd_w = sd_facies(2) * ones(nK, 1);
fprintf('[prior] sd under independence %.4f ; correlation-aware %.4f %.4f %.4f\n', ...
        P.sd_iv, sd_facies(1), sd_facies(2), sd_facies(3));
fprintf('[prior] prior sd used %.4f (rho = %.2f), inflation over IV %.1f%%\n', ...
        sd_facies(2), P.rho_used, 100 * (sd_facies(2) / P.sd_iv - 1));

P.ne_well = ne_w;  P.sd_well = sd_w;  P.sd_facies = sd_facies;
P.xw = xw;  P.yw = yw;  P.Kw = Kw;  P.fw = fw;

% ----------------------------------- correlation length from the measured K --
% The spatial structure of the porosity field is inherited from the facies
% architecture, which the measured log10 K samples directly.  An exponential
% variogram is fitted to the 43 log10 K values by weighted least squares.
lk = log10(Kw);
[hbin, gbin, npair] = fk_empirical_variogram(xw, yw, lk, 14);
ok = npair >= 10 & isfinite(gbin);
sill0 = var(lk);  range0 = 0.3 * max(hbin);
obj = @(p) sum(npair(ok)' .* (gbin(ok)' - ...
          (abs(p(1)) + abs(p(2)) * (1 - exp(-hbin(ok)' / abs(p(3)))))).^2);
p = fminsearch(obj, [0.1*sill0, sill0, range0], ...
               optimset('Display', 'off', 'MaxFunEvals', 5000, 'MaxIter', 5000));
P.nugget = abs(p(1));  P.sill = abs(p(2));  P.vrange = abs(p(3));
P.vario_h = hbin;  P.vario_g = gbin;  P.vario_n = npair;
fprintf(['[prior] exponential variogram on log10 K: nugget %.4f, ' ...
         'partial sill %.4f, range %.0f m\n'], P.nugget, P.sill, P.vrange);

% ------------------------------------------------------------- kriging -----
P.ne = fk_krige(xw, yw, ne_w, G.xa, G.ya, P.nugget, P.sill, P.vrange);
P.sd = fk_krige(xw, yw, sd_w, G.xa, G.ya, P.nugget, P.sill, P.vrange);
P.sd = max(P.sd, min(sd_facies));      % kriging cannot invent a smaller sd

dv = ~isnan(D.DRASTIC);
P.V = fk_krige(D.x(dv), D.y(dv), D.DRASTIC(dv), G.xa, G.ya, ...
               P.nugget, P.sill, P.vrange);
fprintf('[prior] gridded porosity %.4f - %.4f (mean %.4f)\n', ...
        min(P.ne), max(P.ne), mean(P.ne));
fprintf('[prior] gridded DRASTIC  %.1f - %.1f (mean %.1f), from %d wells\n', ...
        min(P.V), max(P.V), mean(P.V), sum(dv));

% ----------------------------------------------------- KL basis of the prior --
% B(x,x') = sd(x) sd(x') exp(-|x-x'|/range).  Nystrom on a strided subset of
% the active cells, then extended to every active cell.
na = G.n_active;
land = 1:cfg.nystrom_stride:na;
nl = numel(land);
Xl = G.xa(land);  Yl = G.ya(land);  Sl = P.sd(land);
Dll = hypot(Xl - Xl', Yl - Yl');
Cll = (Sl .* Sl') .* exp(-Dll / P.vrange);
Cll = (Cll + Cll') / 2;
[U, Lam] = eig(Cll, 'vector');
[Lam, ord] = sort(Lam, 'descend');
U = U(:, ord);
m = min(cfg.n_kl, sum(Lam > 0));
Lam = Lam(1:m);  U = U(:, 1:m);

% Nystrom extension: phi(x) = (1/lambda) * sum_j B(x, x_j) u_j * (nl/na scaling
% is absorbed by renormalising each mode to unit energy in the discrete metric)
Phi = zeros(na, m);
blk = 20000;
for a = 1:blk:na
    b = min(a + blk - 1, na);
    Dxl = hypot(G.xa(a:b) - Xl', G.ya(a:b) - Yl');
    Cxl = (P.sd(a:b) .* Sl') .* exp(-Dxl / P.vrange);
    Phi(a:b, :) = Cxl * (U ./ Lam');
end
for j = 1:m                              % scale so that var captured = Lam
    nrm = sqrt(Phi(:, j)' * Phi(:, j));
    if nrm > 0, Phi(:, j) = Phi(:, j) / nrm; end
end
P.Phi = Phi;
P.lam = Lam * (na / nl);                 % discrete-metric rescaling
P.kl_var_captured = sum(P.lam) / sum(P.sd.^2);
fprintf('[prior] KL basis: %d modes, %.1f%% of the prior variance\n', ...
        m, 100 * P.kl_var_captured);
end

% =========================================================================
function [hb, gb, nb] = fk_empirical_variogram(x, y, z, nbin)
n = numel(x);
[I, J] = find(triu(ones(n), 1));
h = hypot(x(I) - x(J), y(I) - y(J));
g = 0.5 * (z(I) - z(J)).^2;
edges = linspace(0, max(h) * 0.6, nbin + 1);
hb = zeros(nbin, 1);  gb = nan(nbin, 1);  nb = zeros(nbin, 1);
for k = 1:nbin
    s = h >= edges(k) & h < edges(k + 1);
    nb(k) = sum(s);
    if nb(k) > 0
        hb(k) = mean(h(s));
        gb(k) = mean(g(s));
    else
        hb(k) = 0.5 * (edges(k) + edges(k + 1));
    end
end
end

% =========================================================================
function zg = fk_krige(x, y, z, xg, yg, nugget, sill, range)
%FK_KRIGE  Ordinary kriging with an exponential variogram.
n = numel(x);
Dd = hypot(x - x', y - y');
Gm = nugget + sill * (1 - exp(-Dd / range));
Gm(1:n+1:end) = 0;
A = [Gm, ones(n, 1); ones(1, n), 0];
A = A + 1e-10 * eye(n + 1);
zg = zeros(numel(xg), 1);
blk = 5000;
for a = 1:blk:numel(xg)
    b = min(a + blk - 1, numel(xg));
    Dg = hypot(xg(a:b) - x', yg(a:b) - y');
    Gg = nugget + sill * (1 - exp(-Dg / range));
    rhs = [Gg'; ones(1, b - a + 1)];
    w = A \ rhs;
    zg(a:b) = (w(1:n, :)' * z);
end
end
