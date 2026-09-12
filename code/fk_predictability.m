function Pr = fk_predictability(cfg, D)
%FK_PREDICTABILITY  Is the measured nitrate field predictable at all?
%
%  Before any transport model is blamed for a poor cross-validated fit, the
%  question has to be asked of the data.  This routine takes the nitrate
%  values at the wells inside the aquifer and asks how well each of them can
%  be predicted from the others by methods that make no physical assumption
%  whatsoever: inverse-distance weighting, ordinary kriging over a wide range
%  of correlation lengths, and linear regression on the available covariates.
%
%  The number that matters is the leave-one-out coefficient of determination.
%  A value at or below zero means the method does worse than predicting the
%  network mean, which is the honest baseline.

ob = D.inside & D.has_NO3;
x = D.x(ob);  y = D.y(ob);  z = D.NO3(ob);
V = D.DRASTIC(ob);  K = D.K(ob);  dw = D.depth_w(ob);
n = numel(z);
sst = sum((z - mean(z)) .^ 2);
r2 = @(p) 1 - sum((z - p) .^ 2) / sst;
Pr.n = n;

% ------------------------------------------------ inverse distance weights --
Pr.idw_power = [1 2 3];
Pr.idw_R2 = zeros(1, 3);
for q = 1:3
    pw = Pr.idw_power(q);
    pred = zeros(n, 1);
    for i = 1:n
        m = true(n, 1);  m(i) = false;
        d = max(hypot(x(m) - x(i), y(m) - y(i)), 1);
        w = d .^ (-pw);
        pred(i) = sum(w .* z(m)) / sum(w);
    end
    Pr.idw_R2(q) = r2(pred);
    fprintf('[pred] leave-one-out inverse distance, power %d : R2 = %+.3f\n', ...
            pw, Pr.idw_R2(q));
end

% --------------------------------------------------------- ordinary kriging
Pr.krig_range = [1000 2000 4000 8000 16000];
Pr.krig_R2 = zeros(size(Pr.krig_range));
nug = 0.2 * var(z);  sill = 0.8 * var(z);
for q = 1:numel(Pr.krig_range)
    rg = Pr.krig_range(q);
    pred = zeros(n, 1);
    for i = 1:n
        m = true(n, 1);  m(i) = false;
        xs = x(m);  ys = y(m);  zs = z(m);  k = numel(xs);
        Dd = hypot(xs - xs', ys - ys');
        Gm = nug + sill * (1 - exp(-Dd / rg));
        Gm(1:k+1:end) = 0;
        A = [Gm, ones(k, 1); ones(1, k), 0] + 1e-8 * eye(k + 1);
        g0 = nug + sill * (1 - exp(-hypot(xs - x(i), ys - y(i)) / rg));
        w = A \ [g0; 1];
        pred(i) = w(1:k)' * zs;
    end
    Pr.krig_R2(q) = r2(pred);
    fprintf('[pred] leave-one-out ordinary kriging, range %6d m : R2 = %+.3f\n', ...
            rg, Pr.krig_R2(q));
end

% ------------------------------------------------------ linear regressions --
sets = {{'X,Y', [x y]}, ...
        {'X,Y,DRASTIC', [x y V]}, ...
        {'X,Y,DRASTIC,K,depth', [x y V K dw]}, ...
        {'DRASTIC', V}, ...
        {'depth to water', dw}};
Pr.lin_name = cell(1, numel(sets));
Pr.lin_R2 = zeros(1, numel(sets));
for q = 1:numel(sets)
    A = [ones(n, 1), sets{q}{2}];
    pred = zeros(n, 1);
    for i = 1:n
        m = true(n, 1);  m(i) = false;
        b = A(m, :) \ z(m);
        pred(i) = A(i, :) * b;
    end
    Pr.lin_name{q} = sets{q}{1};
    Pr.lin_R2(q) = r2(pred);
    fprintf('[pred] leave-one-out linear on %-22s : R2 = %+.3f\n', ...
            sets{q}{1}, Pr.lin_R2(q));
end

% --------------------------------------------- what the field is organised by
% Reported because the choice of loading model in Section 4.3 turns on it: the
% nitrate pattern follows position far more closely than it follows the
% vulnerability index.
xi_w = x * cosd(90 - 225) + y * sind(90 - 225);
Pr.corr_northing = corr(z, y);
Pr.corr_easting  = corr(z, x);
Pr.corr_alongflow = corr(z, xi_w);
Pr.corr_drastic  = corr(z, V);
Pr.corr_K        = corr(z, K);
Pr.corr_depth    = corr(z, dw);
fprintf(['[pred] nitrate correlates with northing %+.2f, easting %+.2f, ' ...
         'along-flow %+.2f\n'], Pr.corr_northing, Pr.corr_easting, ...
        Pr.corr_alongflow);
fprintf(['[pred] nitrate correlates with DRASTIC %+.2f, K %+.2f, ' ...
         'depth to water %+.2f\n'], Pr.corr_drastic, Pr.corr_K, Pr.corr_depth);

Pr.best_R2 = max([Pr.idw_R2, Pr.krig_R2, Pr.lin_R2]);
fprintf(['[pred] best leave-one-out R2 achieved by any of these methods: ' ...
         '%+.3f\n'], Pr.best_R2);
end
