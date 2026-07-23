clearvars; close all;
%% main.m — Two-region nutrient reaction-diffusion (cell + ECM), spherical
% Diffusivity depends on porosity; porosity is evolved based on bound species.
% Nutrient-ECM Precursor Cycle:
    % a. Coupling cycle
    % b. Diffusion moves nutrients around (depends on current porosity)
    % c. Cells consume nutrients and produce ECM precursors
    % d. ECM precursors link to form bound ECM
    % e. Bound ECM changes the solid volume fraction
    % f. Changed solid volume fraction changes porosity
    % g. Changed porosity changes diffusion coefficients
    % h. Return to step A

%% __ Define Domains __ 
% __ Global __
R_cell = 20;            % cell radius (um)
R_domain = 70;          % domain radius (um)
dr = 1.6;               % spatial step size (um)

% __ Cell Domain: r_start to R_cell __
r_start = 0;  
Nr_cell = round(R_cell/dr);       % number of spatial steps in cell
r_cell_grid = linspace(r_start, R_cell, Nr_cell+1);  % cell radial grid (um)

% __ ECM domain: R_cell to R_domain __
Nr_ecm   = round((R_domain-R_cell)/dr); % number of steps in ECM
r_ecm_grid = linspace(R_cell, R_domain, Nr_ecm+1);  % ECM radial grid (um)

%% __ Parameters __ 
% __ Cell region parameters __ 
k_N_cell = -0.2;       % nutrient reaction rate (1/s) (consumption)
P_N_cell = 0;        % nutrient source term (uM/s) 
k_UM_cell = 0;       % unlinked matrix reaction rate (1/s)
P_UM_cell = 0.1;     % unlinked matrix source term (uM/s)

% __ ECM region parameters __ 
k_N_ecm = 0;        % nutrient reaction rate (1/s)
P_N_ecm = 0;        % nutrient source term (uM/s) 
k_UM_ecm  = -0.048; % unlinked matrix reaction rate (1/s) (consumption due to crosslinking unlinked->linked) (abs(k_m_ecm)=crosslinking rate)
P_UM_ecm  = 0;      % unlinked matrix source term (uM/s)

% __ IC & BC __ 

C0_N_cell = 1;        % initial nutrient concentration (uM), uniform
C0_UM_cell = 0;       % initial matrix concentration (uM), uniform
phi0_S_cell = 0.2;    % initial solid volume fraction
D0_N_cell  = 200;      % initial nutrient diffusivity (um^2/s)
D0_UM_cell  = 1;      % initial unlinked matrix diffusivity (um^2/s)

C0_N_ecm = 1;        % initial nutrient concentration (uM), uniform
C0_UM_ecm = 0;       % initial matrix concentration (uM), uniform
phi0_S_ecm = 0.2;    % initial solid volume fraction
D0_N_ecm  = 2200;      % initial nutrient diffusivity (um^2/s)
D0_UM_ecm  = 1;      % initial unlinked matrix diffusivity (um^2/s)

% __ IC/BC cell __
IC_N_cell   = @(r) C0_N_cell * ones(size(r));   % uniform IC
IC_UM_cell = @(r) C0_UM_cell * ones(size(r));
IC_UM_ecm  = @(r) C0_UM_ecm * ones(size(r));

% __ IC/BC ECM
IC_N_ecm    = @(r) C0_N_ecm * ones(size(r));     % uniform IC
C_N_B  = @(t) C0_N_ecm * ones(size(t));     % far field Boundary: infinite supply
C_UM_B = @(t) C0_UM_ecm * ones(size(t));    % far field Boundary: no matrix supply

% __ General parameters __
t_start = 0;                  % start time (s)
t_end = 200;                  % end time (s)
dt = (2*dr^2) / D0_UM_ecm;    
Nt = round(t_end/dt);         % number of time steps
plot_every_n = max(1, round(Nt / 50));

% __ Define Concentration state arrays __
C_N_cell = zeros(Nr_cell+1, Nt+1);        % cell region solution
C_N_ecm  = zeros(Nr_ecm+1,  Nt+1);        % ECM region solution
C_N_cell(:,1) = IC_N_cell(r_cell_grid)';    % attach ICs
C_N_ecm(:,1)  = IC_N_ecm(r_ecm_grid)';      % attach ICs
C_UM_cell = zeros(Nr_cell+1, Nt+1);  
C_UM_ecm  = zeros(Nr_ecm+1,  Nt+1);  
C_UM_cell(:,1) = IC_UM_cell(r_cell_grid)';
C_UM_ecm(:,1)  = IC_UM_ecm(r_ecm_grid)';

% __ Define Volume Fraction state arrays __
phi_LM_cell = zeros(Nr_cell+1, Nt+1);   % starts at 0 - no crosslinked matrix yet
phi_LM_ecm  = zeros(Nr_ecm+1,  Nt+1);

D_N_cell = D0_N_cell;      
D_UM_cell = D0_UM_cell;      
D_N_ecm = D0_N_ecm;     
D_UM_ecm  = D0_UM_ecm;   
%% __ Time loop  __
for j = 1:Nt
    
  % __ (b) TRANSPORT: diffuse all species  
    % Fick's law continuous flux condition at R_cell
    C_N_interface = interface_cell_ecm(C_N_cell(end-1,j), C_N_ecm(2,j), D_N_cell, D_N_ecm);
    C_UM_interface = interface_cell_ecm(C_UM_cell(end-1,j), C_UM_ecm(2,j), D_UM_cell, D_UM_ecm);

    % Update BCs
    BC_R_n_cell = @(t) C_N_interface * ones(size(t));   % cell surface
    BC_L_n_ecm  = @(t) C_N_interface * ones(size(t));   % ECM inner edge
    BC_L_n_cell = @(t) C_N_cell(2,j) * ones(size(t));   % Zero-flux Neumann BC at cell center (symmetry condition)
    BC_R_m_cell = @(t) C_UM_interface * ones(size(t));
    BC_L_m_ecm  = @(t) C_UM_interface * ones(size(t));
    BC_L_m_cell = @(t) C_UM_cell(2,j) * ones(size(t));

    % Solve Cell region for this time step
    w_temp = crank_nicholson_spherical(r_start, R_cell, ...
        t_start+(j-1)*dt, t_start+j*dt, ...
        Nr_cell, 1, D_N_cell, C_N_cell(2:end-1,j), BC_L_n_cell, ...
        BC_R_n_cell, k_N_cell, P_N_cell, dr, r_cell_grid, ...
        'symmetry');
    C_N_cell(:,j+1) = w_temp(:,end);

    w_temp = crank_nicholson_spherical(r_start, R_cell, ...
        t_start+(j-1)*dt, t_start+j*dt, ...
        Nr_cell, 1, D_UM_cell, C_UM_cell(2:end-1,j), BC_L_m_cell, ...
        BC_R_m_cell, k_UM_cell, P_UM_cell, dr, r_cell_grid, ...
        'symmetry');
    C_UM_cell(:,j+1) = w_temp(:,end);

    % Solve ECM region for this time step
    w_temp = crank_nicholson_spherical(R_cell, R_domain, ...
        t_start+(j-1)*dt, t_start+j*dt, ...
        Nr_ecm, 1, D_N_ecm, C_N_ecm(2:end-1,j), BC_L_n_ecm, ...
        C_N_B, k_N_ecm, P_N_ecm, dr, r_ecm_grid, ...
        'dirichlet');
    C_N_ecm(:,j+1) = w_temp(:,end);

    w_temp = crank_nicholson_spherical(R_cell, R_domain, ...
        t_start+(j-1)*dt, t_start+j*dt, ...
        Nr_ecm, 1, D_UM_ecm, C_UM_ecm(2:end-1,j), BC_L_m_ecm, ...
        C_UM_B, k_UM_ecm, P_UM_ecm, dr, r_ecm_grid, ...
        'dirichlet');
    C_UM_ecm(:,j+1) = w_temp(:,end);

    % __ (c) REACTION: consumption + production 
    %  currently baked into k_N_cell and P_UM_cell as constants.
    %  becomes nutrient-dependent later.

    % __ (d) LINKING: unlinked matrix crosslinks to linked matrix 
    %  phi_LM_cell(:,j+1) = phi_LM_cell(:,j) + ...
    %  phi_LM_ecm(:,j+1)  = phi_LM_ecm(:,j)  + ...

    % __ (e) SOLID FRACTION: update
    %  phi_S_cell = phi_S0_cell + phi_LM_cell(:,j+1);

    % __ (f) POROSITY: update
    %  phi_W_cell = 1 - phi_S_cell;

    % __ (g) DIFFUSIVITY: update 
    %  D_N_cell = D0_N_cell * f(phi_S_cell);   etc.   

end

%% __ Plotting __

% __ Stability print __
fprintf('C at center (r=dr):  %.4f\n', C_N_cell(1, end));
fprintf('C at cell surface:   %.4f\n', C_N_cell(end, end));
fprintf('C at ECM inner edge: %.4f\n', C_N_ecm(1, end));
fprintf('C at far field:      %.4f\n', C_N_ecm(end, end));

specs(1) = struct('name','N cell',  'D',D_N_cell,  'k',k_N_cell,  'L',R_cell);
specs(2) = struct('name','N ecm',   'D',D_N_ecm,   'k',k_N_ecm,   'L',R_domain-R_cell);
specs(3) = struct('name','UM cell', 'D',D_UM_cell, 'k',k_UM_cell, 'L',R_cell);
specs(4) = struct('name','UM ecm',  'D',D_UM_ecm,  'k',k_UM_ecm,  'L',R_domain-R_cell);
fields(1) = struct('name','C_N_cell',  'C',C_N_cell);
fields(2) = struct('name','C_N_ecm',   'C',C_N_ecm);
fields(3) = struct('name','C_UM_cell', 'C',C_UM_cell);
fields(4) = struct('name','C_UM_ecm',  'C',C_UM_ecm);
print_diagnostics(dr, dt, Nt, t_end, specs, fields);

t_vec_plot = linspace(t_start, t_end, Nt+1);

% nutrient: scale 0..max nutrient value
cmin_N = min([C_N_cell(:); C_N_ecm(:)]);
cmax_N = max([C_N_cell(:); C_N_ecm(:)]); 
plot_region_surface(1, r_cell_grid, t_vec_plot, C_N_cell, R_cell, 'Cell', 'Nutrient', cmin_N, cmax_N);
plot_region_surface(2, r_ecm_grid,  t_vec_plot, C_N_ecm,  R_cell, 'ECM',  'Nutrient', cmin_N, cmax_N);

% matrix: scale 0..max matrix value
cmin_UM = min([C_UM_cell(:); C_UM_ecm(:)]);
cmax_UM = max([C_UM_cell(:); C_UM_ecm(:)]);
plot_region_surface(3, r_cell_grid, t_vec_plot, C_UM_cell, R_cell, 'Cell', 'Matrix', cmin_UM, cmax_UM);
plot_region_surface(4, r_ecm_grid,  t_vec_plot, C_UM_ecm,  R_cell, 'ECM',  'Matrix', cmin_UM, cmax_UM);

%% __ Circular heatmap __
% __ Animated circular heatmap __

%animate_cross_section(5, C_N_cell, C_N_ecm, r_cell_grid, r_ecm_grid, ...
%    R_cell, t_vec_plot, dt, plot_every_n, 'Nutrient', cmin_N, cmax_N);
%animate_cross_section(6, C_UM_cell, C_UM_ecm, r_cell_grid, r_ecm_grid, ...
%    R_cell, t_vec_plot, dt, plot_every_n, 'Matrix', cmin_UM, cmax_UM);
plot_snapshot_grid(7, C_N_cell, C_N_ecm, r_cell_grid, r_ecm_grid, R_cell, t_vec_plot, [], 'Nutrient', cmin_N, cmax_N);
plot_snapshot_grid(8, C_UM_cell, C_UM_ecm, r_cell_grid, r_ecm_grid, R_cell, t_vec_plot, [], 'Matrix', cmin_UM, cmax_UM);