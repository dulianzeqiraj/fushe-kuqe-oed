function F = fk_forward(cfg, G, P)
%FK_FORWARD  Transport operator for the Fushe-Kuqe nitrate field.
%
%  Physical setting.  The aquifer is confined and multilayer, recharged from
%  the Mat and Droja rivers and discharging to the Adriatic, with regional flow
%  from northeast to southwest (Cenameri & Beqiraj 2016).  At the published
%  mean seepage velocity the traverse of the aquifer takes of order half a
%  century, which is the same order as the history of intensive agricultural
%  nitrogen loading in the catchment.  The nitrate field is therefore treated
%  as transient rather than steady: it is the response to an areal loading that
%  switched on at an unknown time T0 before the sampling campaign.
%
%  This matters for the whole paper.  At steady state with a uniform Darcy
%  flux the concentration field is almost independent of effective porosity,
%  because n_e cancels between the seepage velocity and the dispersion
%  coefficient.  In the transient problem n_e survives in the storage term and
%  sets how far the front has travelled, so a concentration measurement does
%  carry information about n_e.  The depth-averaged equation solved here is
%
%     n_e dC/dt + q0 dC/dxi = d/dxi (alphaL q0 Gam dC/dxi)
%                           + d/deta (alphaT q0 Gam dC/deta) + m_s(x)
%
%  with q0 the uniform regional Darcy flux, Gam the DRASTIC amplification of
%  dispersion published in Zeqiraj et al. (2026), and m_s the areal loading.
%
%  F is a struct of everything that does not depend on n_e, plus function
%  handles that build the n_e-dependent operator and run the time stepping.

na = G.n_active;

% ------------------------------------------------ uniform regional flux ----
% The workbook records depth to water, not head, and carries no DEM, so the
% hydraulic gradient is not measurable from it.  It is derived from the
% published mean seepage velocity and the measured mean K.
Kmean = mean(P.Kw);
F.i_gradient = cfg.v_mean_pub * P.ne_network_mean / Kmean;
F.q0 = Kmean * F.i_gradient;                       % m/day, uniform Darcy flux
fprintf(['[fwd] mean measured K %.1f m/day, published mean seepage velocity ' ...
         '%.1f m/day\n'], Kmean, cfg.v_mean_pub);
fprintf('[fwd] implied regional gradient %.3e, uniform Darcy flux %.4f m/day\n', ...
        F.i_gradient, F.q0);

% ------------------------------------------- DRASTIC dispersion amplifier --
Vn = (P.V - min(P.V)) / (max(P.V) - min(P.V));
F.Gam = 1 + cfg.alpha_V * Vn .^ cfg.gamma_V;
F.Vn = Vn;
fprintf('[fwd] DRASTIC dispersion amplifier Gamma in [%.2f, %.2f]\n', ...
        min(F.Gam), max(F.Gam));

% ----------------------------------------------- face connectivity table ---
% Neighbour in +xi (downstream, next column) and +eta (across flow, next row).
idx = G.idx;
ny = G.ny;  nx = G.nx;
own = (1:na)';
[ia, ja] = deal(G.ia, G.ja);

nb_xi = zeros(na, 1);
ok = ja < nx;
lin = sub2ind([ny nx], ia(ok), ja(ok) + 1);
nb_xi(ok) = idx(lin);

nb_eta = zeros(na, 1);
ok2 = ia < ny;
lin2 = sub2ind([ny nx], ia(ok2) + 1, ja(ok2));
nb_eta(ok2) = idx(lin2);

F.own = own;  F.nb_xi = nb_xi;  F.nb_eta = nb_eta;
F.n_face_xi  = nnz(nb_xi > 0);
F.n_face_eta = nnz(nb_eta > 0);
fprintf('[fwd] interior faces: %d along flow, %d across flow\n', ...
        F.n_face_xi, F.n_face_eta);

% ---------------------------------------------------- dispersion at faces --
aL = cfg.alphaL;  aT = cfg.alphaT_ratio * cfg.alphaL;
% n_e * D = alpha * q0 * Gamma + n_e * Dm ; the molecular term is six orders
% of magnitude smaller than the mechanical one here and is retained only for
% completeness.
gam_f_xi  = 0.5 * (F.Gam(own(nb_xi > 0))  + F.Gam(nb_xi(nb_xi > 0)));
gam_f_eta = 0.5 * (F.Gam(own(nb_eta > 0)) + F.Gam(nb_eta(nb_eta > 0)));
F.TL = aL * F.q0 * gam_f_xi  * cfg.dy / cfg.dx;   % transmissibility, m^2/day
F.TT = aT * F.q0 * gam_f_eta * cfg.dx / cfg.dy;
F.aL = aL;  F.aT = aT;

% advective transmissibility across a downstream face (upwind, flow in +xi)
F.Aadv = F.q0 * cfg.dy;                            % m^2/day

% ---------------------------------------------------------- source shape ---
F.Vratio = P.V / mean(P.V);

% --------------------------------------------------------- observations ----
F.cell_volume = cfg.dx * cfg.dy;                   % per unit thickness

% ------------------------------------------------------- operator builder --
F.build = @(ne) local_build(cfg, G, F, ne);
F.march = @(A, Mass, s, T0, nsteps) local_march(A, Mass, s, T0, nsteps);
end

% =========================================================================
function [A, Mass] = local_build(cfg, G, F, ne)
%LOCAL_BUILD  Sparse operator A such that  Mass dC/dt = -A C + s.
%
%  A collects the upwind advective flux and the two dispersive fluxes.  Rows
%  are cells; the sign convention is that A C is the net outflow of mass.

na = G.n_active;
I = [];  J = [];  S = [];

% --- dispersion along the flow (xi) ---
k = find(F.nb_xi > 0);
o = F.own(k);  n = F.nb_xi(k);  T = F.TL;
I = [I; o; o; n; n];
J = [J; o; n; n; o];
S = [S; T; -T; T; -T];

% --- dispersion across the flow (eta) ---
k = find(F.nb_eta > 0);
o = F.own(k);  n = F.nb_eta(k);  T = F.TT;
I = [I; o; o; n; n];
J = [J; o; n; n; o];
S = [S; T; -T; T; -T];

% --- upwind advection, flow strictly in +xi ---
% Interior face: mass leaves the upstream cell, enters the downstream one.
k = find(F.nb_xi > 0);
o = F.own(k);  n = F.nb_xi(k);  a = F.Aadv;
I = [I; o; n];
J = [J; o; o];
S = [S; a * ones(numel(o), 1); -a * ones(numel(o), 1)];

% --- outflow faces: cells with no downstream neighbour lose mass freely ---
k = find(F.nb_xi == 0);
I = [I; F.own(k)];
J = [J; F.own(k)];
S = [S; a * ones(numel(k), 1)];

% Inflow faces (no upstream neighbour) carry C = 0, so they add nothing.

A = sparse(I, J, S, na, na);
Mass = spdiags(ne * F.cell_volume, 0, na, na);
end

% =========================================================================
function C = local_march(A, Mass, s, T0, nsteps)
%LOCAL_MARCH  Backward-Euler march from C = 0 to time T0.
dt = T0 / nsteps;
Lhs = decomposition(Mass / dt + A, 'lu');
C = zeros(size(A, 1), 1);
Mdt = Mass / dt;
for it = 1:nsteps
    C = Lhs \ (Mdt * C + s);
end
end
