clear all; close all;

%% __ Define global grid __ 
%(r, dr, Nr, r_start, r_end)
R_cell = 10;            % cell radius (um)
R_domain = 30;          % domain radius (um)
dr = 1;                 % spatial step size (um)
r_start = dr;           % start at dr to avoid singularity at r=0
r_end = R_domain;       % end pos (um)
Nr = round(R/dr);       % number of spatial steps
r = linspace(r_start, r_end, Nr+1);  % radial grid (um)

%% __ Nutrient consumption in cell __ 
% __ Initialize parameters __ 
D_n = 1000;             % diffusion coeff (um^2/s)
k_n = -1;               % reaction rate (1/s) (-) for consumption
P_n = 0;                % source term (uM/s) - no external source
dt_n = 0.0055;          % time step (s)
t_start_n = 0;          % start time (s)
t_end_n = 1;            % end time (s)
Nt_n = round(t_end_n/dt_n);     % number of time steps

% __ IC & BC __
c0_n = 1;               % initial concentration (uM), uniform
IC_n   = @(r) c0_n * ones(size(r));   % u(r,0) = c0, uniform initial concentration
BC_L_n = @(t) c0_n * ones(size(t));   %PLACEHOLDER % u(r_start, t) = c0, nutrients held at cell surface
BC_R_n = @(t) c0_n * ones(size(t));   % u(r_end, t)   = c0, far field assumed infinite source

% __ Call Numerical solver: Crank-Nicholson __ 
figure;
C_n = crank_nicholson_spherical(r_start, r_end, t_start, t_end_n, Nr, Nt_n, D_n, IC_n, BC_L_n, BC_R_n, k_n, P_n, dr, r);

%% __ Matrix production by cell __ 

%% __ Plotting __