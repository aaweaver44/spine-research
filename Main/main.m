clear all; close all;
%% main.m — Two-region nutrient reaction-diffusion (cell + ECM), spherical
%
%   Solves nutrient transport across a cell embedded in ECM, as two coupled
%   1D radial regions sharing a Fick's-law flux interface at r = R_cell:
%
%       Cell:  r = 0 .. R_cell     consumes nutrient (k < 0), center at r = 0
%       ECM:   R_cell .. R_domain  inert medium (k = 0), far field fixed supply
%
%   Structure:
%     - grid & parameters       (per region: D, k, P)
%     - initial / boundary conditions
%     - time loop: interface -> cell solve -> ECM solve
%     - plots: cell surface, ECM surface, circular cross-section
%
%   Calls:  interface_cell_ecm.m,  crank_nicholson_spherical.m
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

%% __ Nutrient consumption in cell: Parameters __ 
% __ Cell region parameters __ 
D_n_cell = 200;              % diffusion coeff inside cell (um^2/s)
k_n_cell = -2;            % reaction rate (1/s) (-) for consumption
P_n_cell = 0;               % source term (uM/s) - no source here

% __ ECM region parameters __ 
D_n_ecm = 2200;             % diffusion coeff in ECM (um^2/s)
k_n_ecm = 0;               % no reaction in ECM
P_n_ecm = 0;                % source term (uM/s) - no source here

% __ General nutrient parameters __
t_start_n = 0;          % start time (s)
t_end_n = 10;            % end time (s)
dt_n = dr^2 / D_n_ecm;  % calculate to set sigma ~1
Nt_n = round(t_end_n/dt_n);     % number of time steps

plot_every = round(Nt_n / 50); 

%% __ Nutrient consumption in cell: IC & BC __ 

c0_n = 1;               % initial concentration (uM), uniform

% __ IC/BC cell __
IC_n_cell   = @(r) c0_n * ones(size(r));   % uniform IC
%BC_L_n_cell = ;   % cell center (placeholder, handled later)
%BC_R_n_cell = ;   % cell surface (placeholder, handled later)

% __ IC/BC ECM
IC_n_ecm    = @(r) c0_n * ones(size(r));   % uniform IC
%BC_L_n_ecm  = ;   % ECM inner edge (placeholder, handled later)
BC_R_n_ecm  = @(t) c0_n * ones(size(t));   % far field: infinite supply

% __ Define solution matrices __
C_n_cell = zeros(Nr_cell+1, Nt_n+1);   % cell region solution
C_n_ecm  = zeros(Nr_ecm+1,  Nt_n+1);   % ECM region solution
C_n_cell(:,1) = IC_n_cell(r_cell_grid)';    % attach ICs
C_n_ecm(:,1)  = IC_n_ecm(r_ecm_grid)';      % attach ICs

%% __ Nutrient consumption in cell: Crank-Nicholson Call __ 
% __ Time loop with interface condition __
for j = 1:Nt_n
    
    % Fick's law continuous flux condition at R_cell
    C_interface = interface_cell_ecm(C_n_cell(end-1,j), C_n_ecm(2,j), D_n_cell, D_n_ecm);

    % Update BCs
    BC_R_n_cell = @(t) C_interface * ones(size(t));   % cell surface
    BC_L_n_ecm  = @(t) C_interface * ones(size(t));   % ECM inner edge
    BC_L_n_cell = @(t) C_n_cell(2,j) * ones(size(t)); % Zero-flux Neumann BC at cell center (symmetry condition)

    % Solve Cell region for this time step
    w_temp = crank_nicholson_spherical(r_start, R_cell, ...
        t_start_n+(j-1)*dt_n, t_start_n+j*dt_n, ...
        Nr_cell, 1, D_n_cell, C_n_cell(2:end-1,j), BC_L_n_cell, ...
        BC_R_n_cell, k_n_cell, P_n_cell, dr, r_cell_grid, ...
        'symmetry');
    C_n_cell(:,j+1) = w_temp(:,end);

    % Solve ECM region for this time step
    w_temp = crank_nicholson_spherical(R_cell, R_domain, ...
        t_start_n+(j-1)*dt_n, t_start_n+j*dt_n, ...
        Nr_ecm, 1, D_n_ecm, C_n_ecm(2:end-1,j), BC_L_n_ecm, ...
        BC_R_n_ecm, k_n_ecm, P_n_ecm, dr, r_ecm_grid, ...
        'dirichlet');
    C_n_ecm(:,j+1) = w_temp(:,end);
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

% __ Nutrient Plots __
t_vec_plot = linspace(t_start_n, t_end_n, Nt_n+1);
cmin = 0;  cmax = c0_n;        % uM range for graph scaling

figure(1);
surf(r_cell_grid, t_vec_plot, C_n_cell')
shading interp;  caxis([cmin cmax]);   % graphs share a scale
xlabel('r (\mum)'); ylabel('t (s)'); zlabel('C (\muM)');
title('Nutrient Concentration - Cell Region');
view(45, 25); 
colorbar;
hold on;
xline(R_cell, 'r-', 'LineWidth', 2, 'Label', 'cell boundary')
hold off;

figure(2);
surf(r_ecm_grid, t_vec_plot, C_n_ecm')
shading interp;  caxis([cmin cmax]);   % graphs share a scale
xlabel('r (\mum)'); ylabel('t (s)'); zlabel('C (\muM)');
title('Nutrient Concentration - ECM Region');
view(45, 25);
colorbar;
hold on;
xline(R_cell, 'r-', 'LineWidth', 2, 'Label', 'cell boundary') 
hold off;

%% __ Circular heatmap __
% __ Animated circular heatmap __

% Combine both grids
r_full = [0, r_cell_grid, r_ecm_grid(2:end)];

% Create 2D circular grid
theta  = linspace(0, 2*pi, 200);
[R_plot, Theta] = meshgrid(r_full, theta);
X = R_plot .* cos(Theta);
Y = R_plot .* sin(Theta);

% Cell boundary circle
theta_circle = linspace(0, 2*pi, 200);
x_circle = R_cell * cos(theta_circle);
y_circle = R_cell * sin(theta_circle);

% Loop through time steps
for j = 1:Nt_n+1
    if mod(j-1, plot_every) ~= 0
        continue
    end

    % Combine solutions at time j
    C_full = [C_n_cell(1,j), C_n_cell(:,j)', C_n_ecm(2:end,j)'];
    C_2D   = repmat(C_full, length(theta), 1);

    % Plot
    figure(3)
    pcolor(X, Y, C_2D);
    shading interp;
    colorbar;
    colormap(jet);
    clim([0 1]);
    axis equal;
    xlabel('x (\mum)'); ylabel('y (\mum)');
    title(sprintf('Nutrient Concentration at t = %.3f s', t_vec_plot(j)));
    hold on;
    plot(x_circle, y_circle, 'w-', 'LineWidth', 1);
    text(0, R_cell+5, 'cell boundary', 'Color', 'w', 'HorizontalAlignment', 'center');
    hold off;

    drawnow;
    pause(0.01);
end