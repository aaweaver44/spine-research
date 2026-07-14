%% Spherical Diffusion-Reaction Animation
%  Solves: dc/dt = D/r^2 * d/dr(r^2 * dc/dr) + k*c
%  Expanded: dc/dt = D*(d2c/dr2 + (2/r)*dc/dr) + k*c
%
%  Scenario: Concentration decaying inward from the cell membrane
%            with a growth reaction term (k > 0)
%
%  Method: Finite differences (explicit forward Euler in time,
%          central differences in space)

clear; close all; clc;

%% ==================== PARAMETERS ====================
% You can modify these to explore different behaviors

D = 1.0;          % Diffusion coefficient (um^2/s)
k = 0.5;          % Reaction rate constant (1/s), positive = growth
R = 10.0;         % Cell radius (um)
Nr = 100;         % Number of spatial grid points
T_final = 5.0;    % Total simulation time (s)

% --- Initial condition: concentration starts high at the membrane ---
% c(r, 0) = gaussian-like profile peaked near r = R

%% ==================== GRID SETUP ====================
dr = R / (Nr - 1);                  % Spatial step size
r = linspace(0, R, Nr)';            % Radial grid (column vector)

% Stability criterion for explicit method (CFL condition)
dt = 0.4 * dr^2 / (2 * D);         % Time step (conservative)
Nt = ceil(T_final / dt);            % Number of time steps
dt = T_final / Nt;                  % Adjust dt to hit T_final exactly

fprintf('Grid: %d spatial points, %d time steps\n', Nr, Nt);
fprintf('dr = %.4f um, dt = %.6f s\n', dr, dt);

%% ==================== INITIAL CONDITION ====================
% Concentration peaked at the membrane (r = R), decaying inward
% Using a Gaussian centered near the membrane
sigma = R * 0.15;                    % Width of initial concentration band
c = exp(-((r - R).^2) / (2 * sigma^2));

% Normalize so max concentration = 1
c = c / max(c);

%% ==================== FIGURE SETUP ====================
figure('Color', 'w', 'Position', [100, 100, 900, 600]);

% --- Main plot: c(r,t) ---
ax1 = subplot(2, 1, 1);
h_line = plot(r, c, 'b-', 'LineWidth', 2.5);
hold on;
h_init = plot(r, c, 'k--', 'LineWidth', 1.0);   % Initial condition reference
hold off;
xlabel('Radial distance r (\mum)', 'FontSize', 13);
ylabel('Concentration c(r,t)', 'FontSize', 13);
title('Spherical Diffusion: Membrane Decay Inward + Growth (k > 0)', 'FontSize', 14);
legend('Current c(r,t)', 'Initial condition', 'Location', 'northwest');
xlim([0, R]);
ylim([0, 2.0]);  % Allow room for growth
grid on;
set(gca, 'FontSize', 12);

% --- Time annotation ---
h_time = text(0.5, 1.8, 't = 0.000 s', 'FontSize', 14, 'FontWeight', 'bold');

% --- Lower plot: 2D cross-section visualization ---
ax2 = subplot(2, 1, 2);
theta_vis = linspace(0, 2*pi, 200);
[R_mesh, Theta_mesh] = meshgrid(r, theta_vis);
C_mesh = repmat(c', length(theta_vis), 1);   % Same c at all angles (symmetry)
X_vis = R_mesh .* cos(Theta_mesh);
Y_vis = R_mesh .* sin(Theta_mesh);

h_surf = pcolor(X_vis, Y_vis, C_mesh);
shading interp;
colormap(ax2, hot);
colorbar;
caxis([0, 2.0]);
axis equal;
xlim([-R, R] * 1.1);
ylim([-R, R] * 1.1);
xlabel('x (\mum)', 'FontSize', 13);
ylabel('y (\mum)', 'FontSize', 13);
title('Cross-Section of Cell (Spherical Symmetry)', 'FontSize', 14);
set(gca, 'FontSize', 12);

% Draw cell membrane circle
hold on;
plot(R*cos(theta_vis), R*sin(theta_vis), 'w-', 'LineWidth', 2);
hold off;

%% ==================== TIME-STEPPING LOOP ====================
% Using explicit finite difference method
%
% Interior points (i = 2 to Nr-1):
%   dc/dt = D * [ (c(i+1) - 2c(i) + c(i-1))/dr^2
%               + (2/r(i)) * (c(i+1) - c(i-1))/(2*dr) ]
%           + k * c(i)
%
% Boundary conditions:
%   r = 0:  symmetry  =>  dc/dr = 0  (use L'Hopital: Laplacian = 3D * d2c/dr2)
%   r = R:  fixed concentration at membrane  =>  c(Nr) = exp(-k_decay * t)
%           (membrane concentration slowly depleting as molecules diffuse inward)

% Animation frame rate
frames_per_plot = max(1, floor(Nt / 300));  % ~300 frames total

% Store membrane decay rate (membrane concentration decays over time)
k_membrane = 0.3;  % Rate at which membrane source depletes

for n = 1:Nt
    t_current = n * dt;
    c_new = c;
    
    % --- Interior points: i = 2 to Nr-1 ---
    for i = 2:(Nr-1)
        % Central difference for d2c/dr2
        d2c_dr2 = (c(i+1) - 2*c(i) + c(i-1)) / dr^2;
        
        % Central difference for dc/dr
        dc_dr = (c(i+1) - c(i-1)) / (2 * dr);
        
        % Spherical Laplacian: d2c/dr2 + (2/r)*dc/dr
        laplacian = d2c_dr2 + (2 / r(i)) * dc_dr;
        
        % Full equation: dc/dt = D * laplacian + k * c
        c_new(i) = c(i) + dt * (D * laplacian + k * c(i));
    end
    
    % --- Boundary at r = 0 (symmetry): L'Hopital's rule ---
    % As r -> 0, (2/r)(dc/dr) -> 2*(d2c/dr2), so Laplacian -> 3*d2c/dr2
    d2c_dr2_center = 2 * (c(2) - c(1)) / dr^2;  % Forward difference
    laplacian_center = 3 * d2c_dr2_center;
    c_new(1) = c(1) + dt * (D * laplacian_center + k * c(1));
    
    % --- Boundary at r = R (membrane): decaying source ---
    c_new(Nr) = exp(-k_membrane * t_current);
    
    % Prevent negative concentrations
    c_new = max(c_new, 0);
    
    % Update solution
    c = c_new;
    
    % --- Update animation ---
    if mod(n, frames_per_plot) == 0 || n == Nt
        % Update line plot
        set(h_line, 'YData', c);
        set(h_time, 'String', sprintf('t = %.3f s', t_current));
        
        % Update 2D cross-section
        C_mesh = repmat(c', length(theta_vis), 1);
        set(h_surf, 'CData', C_mesh);
        
        drawnow;
        pause(0.01);  % Small pause for visible animation
    end
end

fprintf('Simulation complete! Final time: %.3f s\n', T_final);
fprintf('Peak concentration: %.4f at r = %.2f um\n', max(c), r(find(c == max(c), 1)));

%% ==================== EXPLANATION ====================
% What you should see:
%
% 1. Initially, concentration is high at the membrane (r = R) and 
%    near zero at the center (r = 0).
%
% 2. Over time, diffusion causes the concentration to spread INWARD
%    toward the center of the cell.
%
% 3. The growth term (k > 0) amplifies the concentration everywhere,
%    competing with the diffusion spreading.
%
% 4. The membrane source slowly depletes, so eventually the inward
%    diffusion + growth dynamics dominate.
%
% The 2/r correction term from spherical coordinates means diffusion
% behaves differently than in a flat geometry -- there is a geometric
% focusing effect as molecules move toward the center (smaller shells).