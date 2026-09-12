function R = oed_local(Phi, lam, w, sd_loc, cand, nsel, xa, ya, min_sep)
if nargin < 9, min_sep = 0; end
%OED_LOCAL  D-optimal and goal-oriented placement for a local observation.
%
%  A characterisation point placed at x measures the field there directly, so
%  its Jacobian row in the KL basis is simply phi(x) and the two criteria take
%  the same rank-one form as for a transport observation:
%
%     D-optimal    log(1 + phi Sigma phi' / sd^2)
%     goal         (c' Sigma phi')^2 / sd^2 / (1 + phi Sigma phi' / sd^2)
%
%  with c the coefficient vector of the vulnerability-weighted mean of the
%  field, g = sum_x w(x) dn(x) / sum_x w(x).
%
%  Using the same operator at both sites is what makes the two comparable:
%  it removes the transport model, and with it every site-specific assumption
%  about flow, from the comparison.

Sigma0 = diag(lam);
c = (Phi' * w) / sum(w);
R.c = c;
R.var_prior = c' * Sigma0 * c;

for mode = {'D', 'goal'}
    md = mode{1};
    Sg = Sigma0;
    cc = cand;
    sel = zeros(nsel, 1);
    vg = zeros(nsel + 1, 1);  vg(1) = R.var_prior;
    ld = zeros(nsel, 1);
    for j = 1:nsel
        SgC = Sg * c;
        PS = Phi * Sg;
        quad = sum(PS .* Phi, 2);
        switch md
            case 'D'
                score = log(1 + quad / sd_loc^2);
            case 'goal'
                score = ((Phi * SgC) .^ 2 / sd_loc^2) ./ (1 + quad / sd_loc^2);
        end
        score(~cc) = -inf;
        [~, b] = max(score);
        sel(j) = b;  cc(b) = false;
        if min_sep > 0
            cc(hypot(xa - xa(b), ya - ya(b)) < min_sep) = false;
        end
        pb = Phi(b, :);
        u = Sg * pb';
        ld(j) = log(1 + (pb * u) / sd_loc^2);
        Sg = Sg - (u * u') / (sd_loc^2 + pb * u);
        Sg = (Sg + Sg') / 2;
        vg(j + 1) = c' * Sg * c;
    end
    R.(md).cells = sel;
    R.(md).var_goal = vg;
    R.(md).logdet_gain = ld;
end

% how far apart the two designs are, and how they differ in vulnerability
d = zeros(nsel, 1);
for j = 1:nsel
    d(j) = min(hypot(xa(R.goal.cells(j)) - xa(R.D.cells), ...
                     ya(R.goal.cells(j)) - ya(R.D.cells)));
end
R.separation_m = d;
R.median_separation_m = median(d);
R.var_reduction_D    = 1 - R.D.var_goal(end)    / R.var_prior;
R.var_reduction_goal = 1 - R.goal.var_goal(end) / R.var_prior;
R.goal_advantage = (R.D.var_goal(end) - R.goal.var_goal(end)) / R.D.var_goal(end);
end
