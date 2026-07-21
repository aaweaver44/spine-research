clearvars; close all;
%% main.m — Two-region nutrient reaction-diffusion (cell + ECM), spherical
%
%   Solves nutrient and matrix transport across a cell embedded in ECM
%       Cell:  r = 0 .. R_cell     
%         - nutrient: consumes (k < 0)
%         - matrix: produces (k > 0)
%       ECM:   R_cell .. R_domain  
%         - nutrient: inert medium (k = 0), far field fixed supply
%
%   Structure:
%     - grid & parameters       (per region: D, k, P)
%     - initial / boundary conditions
%     - time loop: interface -> cell solve -> ECM solve
%     - plots: cell surface, ECM surface, circular cross-section
%
%% __ Define Domains __ 
% __ Global __
R_cell = 20;            % cell radius (um)
R_domain = 70;          % domain radius (um)
dr = 2;                 % spatial step size (um)

% __ Cell Domain: r_start to R_cell __
r_start = 0;  
Nr_cell = round(R_cell/dr);       % number of spatial steps in cell
r_cell_grid = linspace(r_start, R_cell, Nr_cell+1);  % cell radial grid (um)

% __ ECM domain: R_cell to R_domain __
Nr_ecm   = round((R_domain-R_cell)/dr); % number of steps in ECM
r_ecm_grid = linspace(R_cell, R_domain, Nr_ecm+1);  % ECM radial grid (um)

%% __ Parameters __ 
% __ Cell region parameters __ 
D_n_cell = 200;     % nutrient diffusion coeff inside cell (um^2/s)
k_n_cell = -5;      % nutrient reaction rate (1/s) (consumption)
P_n_cell = 0;       % nutrient source term (uM/s) 
D_m_cell = 50;      % matrix diffusivity in cell (um^2/s)
k_m_cell = 0;       % matrix reaction rate (1/s)
P_m_cell = 0.1;    % matrix source term (uM/s)

% __ ECM region parameters __ 
D_n_ecm = 2200;     % nutrient diffusion coeff in ECM (um^2/s)
k_n_ecm = 0;        % nutrient reaction rate (1/s)
P_n_ecm = 0;        % nutrient source term (uM/s) 
D_m_ecm  = 300;     % matrix diffusivity in ECM (um^2/s)
k_m_ecm  = -0.05;   % matrix reaction rate (1/s) (consumption due to crosslinking unbound->bound)
P_m_ecm  = 0;       % matrix source term (uM/s)

% __ General parameters __
t_start_n = 0;          % start time (s)
t_end_n = 10;            % end time (s)
dt_n = dr^2 / D_n_ecm;  % calculate to set sigma ~1
Nt_n = round(t_end_n/dt_n);     % number of time steps

plot_every_n = max(1, round(Nt_n / 50));

%% __ IC & BC __ 

c0_n = 1;       % nutrient initial concentration (uM), uniform
c0_m = 0;       % matrix initial concentration (uM), uniform

% __ IC/BC cell __
IC_n_cell   = @(r) c0_n * ones(size(r));   % uniform IC
IC_m_cell = @(r) c0_m * ones(size(r));
IC_m_ecm  = @(r) c0_m * ones(size(r));

% __ IC/BC ECM
IC_n_ecm    = @(r) c0_n * ones(size(r));   % uniform IC
BC_R_n_ecm  = @(t) c0_n * ones(size(t));   % far field: infinite supply
BC_R_m_ecm = @(t) c0_m * ones(size(t));    % far field: no matrix supply

% __ Define solution matrices __
C_n_cell = zeros(Nr_cell+1, Nt_n+1);   % cell region solution
C_n_ecm  = zeros(Nr_ecm+1,  Nt_n+1);   % ECM region solution
C_n_cell(:,1) = IC_n_cell(r_cell_grid)';    % attach ICs
C_n_ecm(:,1)  = IC_n_ecm(r_ecm_grid)';      % attach ICs
C_m_cell = zeros(Nr_cell+1, Nt_n+1);  
C_m_ecm  = zeros(Nr_ecm+1,  Nt_n+1);  
C_m_cell(:,1) = IC_m_cell(r_cell_grid)';
C_m_ecm(:,1)  = IC_m_ecm(r_ecm_grid)';

%% __ Crank-Nicholson Call __ 
% __ Time loop with interface condition __
for j = 1:Nt_n
    
    % Fick's law continuous flux condition at R_cell
    C_n_interface = interface_cell_ecm(C_n_cell(end-1,j), C_n_ecm(2,j), D_n_cell, D_n_ecm);
    C_m_interface = interface_cell_ecm(C_m_cell(end-1,j), C_m_ecm(2,j), D_m_cell, D_m_ecm);

    % Update BCs
    BC_R_n_cell = @(t) C_n_interface * ones(size(t));   % cell surface
    BC_L_n_ecm  = @(t) C_n_interface * ones(size(t));   % ECM inner edge
    BC_L_n_cell = @(t) C_n_cell(2,j) * ones(size(t));   % Zero-flux Neumann BC at cell center (symmetry condition)
    BC_R_m_cell = @(t) C_m_interface * ones(size(t));
    BC_L_m_ecm  = @(t) C_m_interface * ones(size(t));
    BC_L_m_cell = @(t) C_m_cell(2,j) * ones(size(t));

    % Solve Cell region for this time step
    w_temp = crank_nicholson_spherical(r_start, R_cell, ...
        t_start_n+(j-1)*dt_n, t_start_n+j*dt_n, ...
        Nr_cell, 1, D_n_cell, C_n_cell(2:end-1,j), BC_L_n_cell, ...
        BC_R_n_cell, k_n_cell, P_n_cell, dr, r_cell_grid, ...
        'symmetry');
    C_n_cell(:,j+1) = w_temp(:,end);

    w_temp = crank_nicholson_spherical(r_start, R_cell, ...
        t_start_n+(j-1)*dt_n, t_start_n+j*dt_n, ...
        Nr_cell, 1, D_m_cell, C_m_cell(2:end-1,j), BC_L_m_cell, ...
        BC_R_m_cell, k_m_cell, P_m_cell, dr, r_cell_grid, ...
        'symmetry');
    C_m_cell(:,j+1) = w_temp(:,end);

    % Solve ECM region for this time step
    w_temp = crank_nicholson_spherical(R_cell, R_domain, ...
        t_start_n+(j-1)*dt_n, t_start_n+j*dt_n, ...
        Nr_ecm, 1, D_n_ecm, C_n_ecm(2:end-1,j), BC_L_n_ecm, ...
        BC_R_n_ecm, k_n_ecm, P_n_ecm, dr, r_ecm_grid, ...
        'dirichlet');
    C_n_ecm(:,j+1) = w_temp(:,end);

    w_temp = crank_nicholson_spherical(R_cell, R_domain, ...
        t_start_n+(j-1)*dt_n, t_start_n+j*dt_n, ...
        Nr_ecm, 1, D_m_ecm, C_m_ecm(2:end-1,j), BC_L_m_ecm, ...
        BC_R_m_ecm, k_m_ecm, P_m_ecm, dr, r_ecm_grid, ...
        'dirichlet');
    C_m_ecm(:,j+1) = w_temp(:,end);

end

fprintf('C at center (r=dr):  %.4f\n', C_n_cell(1, end));
fprintf('C at cell surface:   %.4f\n', C_n_cell(end, end));
fprintf('C at ECM inner edge: %.4f\n', C_n_ecm(1, end));
fprintf('C at far field:      %.4f\n', C_n_ecm(end, end));

%% __ Plotting __

% __ Stability print __
sigma_n_cell = D_n_cell * dt_n / dr^2;
sigma_n_ecm  = D_n_ecm  * dt_n / dr^2;
rho_n_cell   = k_n_cell * dt_n;
fprintf('\n-- Nutrients --\n')
fprintf('sigma cell = %f\n', sigma_n_cell);
fprintf('sigma ECM  = %f\n', sigma_n_ecm);
fprintf('rho cell   = %f\n', rho_n_cell);

t_vec_plot = linspace(t_start_n, t_end_n, Nt_n+1);

% nutrient: scale 0..max nutrient value
cmax_n = max([C_n_cell(:); C_n_ecm(:)]); 
plot_region_surface(1, r_cell_grid, t_vec_plot, C_n_cell, R_cell, 'Cell', 'Nutrient', cmax_n);
plot_region_surface(2, r_ecm_grid,  t_vec_plot, C_n_ecm,  R_cell, 'ECM',  'Nutrient', cmax_n);

% matrix: scale 0..max matrix value
cmax_m = max([C_m_cell(:); C_m_ecm(:)]);
plot_region_surface(3, r_cell_grid, t_vec_plot, C_m_cell, R_cell, 'Cell', 'Matrix', cmax_m);
plot_region_surface(4, r_ecm_grid,  t_vec_plot, C_m_ecm,  R_cell, 'ECM',  'Matrix', cmax_m);

%% __ Circular heatmap __
% __ Animated circular heatmap __

cmax_n = max([C_n_cell(:); C_n_ecm(:)]);
cmax_m = max([C_m_cell(:); C_m_ecm(:)]);
animate_cross_section(5, C_n_cell, C_n_ecm, r_cell_grid, r_ecm_grid, ...
    R_cell, t_vec_plot, dt_n, plot_every_n, 'Nutrient', cmax_n);
animate_cross_section(6, C_m_cell, C_m_ecm, r_cell_grid, r_ecm_grid, ...
    R_cell, t_vec_plot, dt_n, plot_every_n, 'Matrix', cmax_m);
plot_snapshot_grid(7, C_n_cell, C_n_ecm, r_cell_grid, r_ecm_grid, R_cell, t_vec_plot, [], 'Nutrient', cmax_n);
plot_snapshot_grid(8, C_m_cell, C_m_ecm, r_cell_grid, r_ecm_grid, R_cell, t_vec_plot, [], 'Matrix',   cmax_m);