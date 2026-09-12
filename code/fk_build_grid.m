function G = fk_build_grid(cfg, D)
%FK_BUILD_GRID  Flow-aligned regular grid over the aquifer polygon.
%
%  The grid axes are rotated onto the regional flow direction reported by
%  Cenameri & Beqiraj (2016): northeast to southwest.  Aligning the grid with
%  the flow makes the dispersion tensor diagonal in grid coordinates, which
%  removes the cross-derivative terms from the transport operator and with
%  them the main source of discretisation error in a rotated-flow problem.
%
%  Fields
%    G.exi, G.eta   unit vectors of the rotated frame in real-world coordinates
%    G.xa,  G.ya    real-world coordinates of the active cell centres
%    G.ia,  G.ja    rotated-frame (row, column) index of each active cell
%    G.idx          (ny x nx) map from grid position to active-cell number

az = cfg.azim_flow_deg;
th = deg2rad(90 - az);                 % compass azimuth to mathematical angle
exi = [cos(th), sin(th)];              % unit vector along the flow
eta = [-sin(th), cos(th)];             % unit vector across the flow

% rotated coordinates of the boundary polygon
xi_b  = D.bx * exi(1) + D.by * exi(2);
eta_b = D.bx * eta(1) + D.by * eta(2);

xi0 = min(xi_b);  xi1 = max(xi_b);
et0 = min(eta_b); et1 = max(eta_b);

nx = ceil((xi1 - xi0) / cfg.dx);       % columns advance along the flow
ny = ceil((et1 - et0) / cfg.dy);       % rows advance across the flow

xic = xi0 + (0.5:1:(nx - 0.5)) * cfg.dx;
etc = et0 + (0.5:1:(ny - 0.5)) * cfg.dy;
[XI, ET] = meshgrid(xic, etc);

% back to real-world coordinates (the frame is orthonormal, so the inverse is
% the transpose)
XX = XI * exi(1) + ET * eta(1);
YY = XI * exi(2) + ET * eta(2);

in = inpolygon(XX, YY, D.bx, D.by);

G = struct();
G.exi = exi;  G.eta = eta;  G.theta = th;
G.nx = nx;  G.ny = ny;
G.xi = xic;  G.eta_c = etc;
G.X = XX;   G.Y = YY;
G.mask = in;
G.idx = zeros(ny, nx);
G.idx(in) = 1:nnz(in);
G.n_active = nnz(in);
G.cell_area = cfg.dx * cfg.dy;
G.area_km2 = G.n_active * G.cell_area / 1e6;
G.xa = XX(in);
G.ya = YY(in);
[G.ia, G.ja] = find(in);               % row (across flow), column (along flow)

% area of the digitised polygon itself, for the discretisation check
G.area_polygon_km2 = polyarea(D.bx, D.by) / 1e6;

fprintf(['[grid] flow azimuth %g deg, rotated grid %d x %d cells of %g m\n'], ...
        az, ny, nx, cfg.dx);
fprintf(['[grid] %d active cells, discretised area %.2f km2 ' ...
         '(polygon %.2f km2, difference %.2f%%)\n'], ...
        G.n_active, G.area_km2, G.area_polygon_km2, ...
        100 * abs(G.area_km2 - G.area_polygon_km2) / G.area_polygon_km2);
end
