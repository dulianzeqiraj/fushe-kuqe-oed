function cfg = fk_config()
%FK_CONFIG  Every constant used by the pipeline, with its provenance.
%
%  Nothing in this file is tuned to reproduce a target result.  Each value is
%  either (a) measured and present in data/, (b) taken from a cited published
%  source, or (c) declared here as a modelling assumption and swept in the
%  sensitivity analysis (see run_all.m, stage S).
%
%  Provenance tags:
%    [DATA]  read from the field workbook, nothing assumed
%    [PUB]   published value, citation given
%    [ASSUM] modelling assumption made in this paper, swept in sensitivity
%    [DERIV] derived here from [DATA] and [PUB] quantities, formula given

cfg = struct();

% ---------------------------------------------------------------- paths ---
here = fileparts(mfilename('fullpath'));
cfg.dir_code    = here;
cfg.dir_data    = fullfile(here, '..', 'data');
cfg.dir_results = fullfile(here, '..', 'results');
cfg.dir_figures = fullfile(here, '..', 'figures');

% ------------------------------------------------------------ numerical ---
cfg.dx = 100;                 % [ASSUM] cell size, m.  Advection is first-order
cfg.dy = 100;                 %         upwind, whose numerical dispersion is
                              %         equivalent to a dispersivity of dx/2 =
                              %         50 m, ten per cent of the baseline
                              %         alphaL of 500 m.  The dispersivity
                              %         sweep spans that error several times
                              %         over, so the conclusions do not rest
                              %         on the discretisation.
cfg.n_kl = 150;              % [ASSUM] Karhunen-Loeve modes retained
cfg.nystrom_stride = 8;       % [ASSUM] coarse-subgrid stride for Nystrom

% ------------------------------------------------ aquifer / flow physics ---
cfg.b_aquifer = 40;           % [PUB] depth-averaged saturated thickness, m.
                              %       Zeqiraj et al. 2026, J Hazard Mater Adv
                              %       23, 101261.  The primary hydrogeological
                              %       description (Cenameri & Beqiraj 2016)
                              %       gives a confined multilayer system whose
                              %       gravel-sand thickness grows from 5-10 m
                              %       in the east to 180-200 m in the west; the
                              %       single value is a depth-averaged 2-D
                              %       idealisation and is stated as such.
cfg.v_mean_pub = 2.1;         % [PUB] mean seepage velocity, m/day, same source
cfg.alphaL = 500;             % [ASSUM] longitudinal dispersivity, m, from the
                              %         Gelhar et al. (1992) scale dependence
                              %         alphaL ~ 0.1 L at a transport scale of
                              %         a few km.  Swept over 100-1000 m.
cfg.alphaT_ratio = 0.1;       % [PUB] alphaT = alphaL/10, Zheng & Bennett 2002
cfg.Dm = 1.9e-9 * 86400;      % [PUB] molecular diffusion of nitrate, m^2/day
                              %       (1.9e-9 m^2/s, Li & Gregory 1974)

% Regional flow direction.  The primary hydrogeological description of this
% aquifer (Cenameri & Beqiraj 2016, Bull. Geol. Soc. Greece 50(2), 665,
% doi:10.12681/bgsg.11772) gives flow from the Mat river outlet in the
% northeast to the Adriatic discharge boundary in the southwest, with the
% seawater wedge advancing in the opposite sense.  A single regional azimuth
% is therefore used rather than the rotating field of earlier work.
cfg.azim_flow_deg = 225;      % [PUB] flow towards SW (NE -> SW)
% Gradient magnitude is not measured: the workbook records depth to water, not
% head, and no DEM accompanies it.  It is derived from the published mean
% seepage velocity and the measured mean K as i = v n_e / K.          [DERIV]
cfg.i_gradient = [];          % filled in fk_flow_field from measured K

% ------------------------------------------------- DRASTIC-transport link --
% Vulnerability amplification of the dispersion tensor, the form published in
% Zeqiraj et al. 2026 (J Hazard Mater Adv 23, 101261), Eq. (6).
cfg.alpha_V = 0.5;            % [ASSUM] amplitude, swept 0 - 1
cfg.gamma_V = 1.0;            % [ASSUM] exponent, swept 0.5 - 2

% ------------------------------------------------------- porosity prior ----
% Facies reference ranges for effective porosity, and the pathway error
% correlations, are taken from the fusion already published for this aquifer:
%   Zeqiraj & Beqiraj 2026, J Contam Hydrol 283, 105086 (Table 1, Eq. 7)
% This paper does not re-derive them; it consumes them as the prior.   [PUB]
cfg.facies_names = {'Gravel - coarse-grained', ...
                    'Gravel - medium-grained', ...
                    'Gravel - fine-grained'};
cfg.facies_short = {'coarse', 'medium', 'fine'};
cfg.ne_lo   = [0.30, 0.28, 0.25];   % [PUB] lower end of reference range
cfg.ne_hi   = [0.34, 0.32, 0.29];   % [PUB] upper end of reference range
cfg.rho_pathway = [0.15, 0.35, 0.58];  % [PUB] facies pathway-error correlation
cfg.sigma_VS = 0.042;         % [PUB] Vukovic-Soro pathway sd
cfg.sigma_KC = 0.051;         % [PUB] Kozeny-Carman pathway sd

% ------------------------------------------------------- observations -----
cfg.sigma_obs = 1.0;          % [ASSUM] nitrate measurement sd, mg/L.  Swept
                              %         0.5 - 2.0 mg/L.
cfg.sigma_K_rel = 0.10;       % [ASSUM] relative sd of a pumping-test K, used
                              %         only for the information-equivalence
                              %         comparison.  Swept 0.05 - 0.20.

% -------------------------------------------------------- calibration -----
cfg.cal_dt = 200;             % time step of the calibration march, days
cfg.cal_nsteps = 400;         % 80 000 days = 219 yr horizon, comfortably
                              % longer than any plausible loading history
cfg.cal_kappas = [0 0.25 0.5 0.75 1 1.5 2 3];   % DRASTIC loading exponents

% ------------------------------------------------------------- design -----
cfg.n_new_wells = 5;          % how many augmentation sites to place
cfg.gamma_Vdesign = 1.5;      % [ASSUM] vulnerability weighting exponent in
                              %         the D_V criterion.  Swept 0 - 3.
cfg.horizon_frac = 0.90;      % fraction of sensitivity-kernel mass defining
                              %         the information horizon
cfg.sd_local = 0.02;          % [ASSUM] sd of a direct porosity determination
                              %         at a characterisation point, used by
                              %         the local observation operator that
                              %         both sites share
cfg.mrva_dx = 1000;           % cell size of the MRVA grid, m (the well
                              %         spacing there is tens of km)

% -------------------------------------------------------- particle step ---
cfg.gn_iters = 12;            % Gauss-Newton iterations for the MAP estimate
cfg.max_desc = 20000;         % cap on the descent-comparison iteration count
cfg.M_particles = 50;         % primary implicit-sampling particles
cfg.M_sir = 500;              % particles for the SIR comparison
cfg.seed = 20260912;          % fixed seed, reported in the paper

end
