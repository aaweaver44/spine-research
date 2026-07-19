clear all; close all;
%   Solves: dc/dt = D/r^2 * d/dr(r^2 * dc/dr) + k*c + rho
%           = D * [d2c/dr2 + (2/r)*dc/dr] + k*c + rho
%
%   Terms:
%      Domain: r = 0 (center of cell) to r = R (cell membrane)
%      D * [...]  = diffusion (spatial spreading through the sphere)
%      k * c      = reaction  (growth/decay proportional to current concentration)
%      rho        = constant source (fixed production rate everywhere, independent of c)

%% %__Initialize parameters, grid, & ICs for spherical diffusion__%%%

D = 1;         % diffusion coefficient (cm^2/s)
k = 1;         % reaction rate (1/s), positive = growth
rho = 0.5;     % constant source term (concentration/s)
               %   positive = production, negative = consumption
               %   try rho = 0 to recover the original equation without it
R = 1;         % cell radius (cm)
dr = 0.1;      % spatial step size (cm)

dt = 0.003;    % time step size (sec)
               %   good value for crank:   0.0055
               %   good value for forward: 0.003

t_end = 1;     % end simulation time (sec)

r_start = 0;             % center of cell
r_end   = r_start + R;   % membrane

Nr = round(R / dr);                          % number of spatial steps

r = linspace(r_start, r_end, Nr+1);          % spatial grid
t_start = 0;
Nt = round(t_end / dt);                      % number of time steps
t_matrix = linspace(t_start, t_end, Nt+1);   % time grid

%%%__Boundary & Initial Conditions__%%%
c0  = 1;     % initial concentration scale
rbc = 0;     % membrane boundary value

% Initial condition: sine-squared profile
IC = @(r) sin(2*pi*r/R).^2;

% Boundary at r = R (membrane): constant value
BC_R = @(t) rbc * ones(size(t));

%% %__Numerical Solutions__%%%

% Crank-Nicolson (implicit, spherical)
figure;
crank_matrix = crank_nicholson_spherical_p(r_start, r_end, t_start, t_end, ...
                                          Nr, Nt, D, IC, BC_R, k, rho);

%% %__Forward Euler Comparison__%%%
% Forward Euler (explicit, spherical)
figure;
forward_matrix = forward_difference_spherical_p(r_start, r_end, t_start, t_end, ...
                                              Nr, Nt, D, IC, BC_R, k, rho);
%%%__Compare the two methods at final time (line plot)__%%%
figure;
r_full = linspace(r_start, r_end, Nr+1);
plot(r_full, crank_matrix(:,end), 'b-o', 'LineWidth', 1.5, 'DisplayName', 'Crank-Nicolson');
hold on;
plot(r_full, forward_matrix(:,end), 'r--x', 'LineWidth', 1.5, 'DisplayName', 'Forward Difference');
xlabel('r (radius)');
ylabel('Concentration c(r, t_{end})');
title(sprintf('Comparison at Final Time  (\\rho = %.2f)', rho));
legend('Location', 'best');
grid on;
hold off;

%% %==========================================================%%%
%%%                     Visualization                        %%%
%%%==========================================================%%%

% Select solution to animate
sol = crank_matrix;   
method_name = 'Crank-Nicolson';

% Build polar mesh for the circular cross-section
N_theta = 200;
theta = linspace(0, 2*pi, N_theta);
[R_mesh, Theta_mesh] = meshgrid(r_full, theta);
X = R_mesh .* cos(Theta_mesh);
Y = R_mesh .* sin(Theta_mesh);

% Color limits
c_min = min(sol(:));
c_max = max(sol(:));
if c_min == c_max
    c_max = c_min + 1;
end

% ---- Figure: Snapshot grid (6 time slices) ----
figure('Color', 'w', 'Position', [50, 50, 1200, 700]);
sgtitle(sprintf('Cell Cross-Section Over Time  (%s,  \\rho = %.2f)', method_name, rho), ...
        'FontSize', 16, 'FontWeight', 'bold');

num_snapshots = 6;
% linear spacing
% snapshot_indices = round(linspace(1, size(sol, 2), num_snapshots));
% logarithmic spacing
snapshot_indices = unique(round(logspace(0, log10(size(sol,2)), num_snapshots)));

for s = 1:num_snapshots
    j = snapshot_indices(s);
    t_now = (j-1) * dt;
    
    c_snapshot = sol(:, j)';
    C_mesh = repmat(c_snapshot, N_theta, 1);
    
    subplot(2, 3, s);
    pcolor(X, Y, C_mesh);
    shading interp;
    colormap(hot);
    caxis([c_min, c_max]);
    axis equal tight off;
    
    hold on;
    plot(R*cos(theta), R*sin(theta), 'w-', 'LineWidth', 1.5);
    hold off;
    
    title(sprintf('t = %.3f s', t_now), 'FontSize', 12);
end

cb = colorbar('Position', [0.93, 0.15, 0.02, 0.7]);
cb.Label.String = 'Concentration';
cb.Label.FontSize = 12;

% ---- Figure: Full animation ----
fig_anim = figure('Color', 'w', 'Position', [100, 100, 800, 700]);

c_init = sol(:, 1)';
C_mesh = repmat(c_init, N_theta, 1);

h_pcolor = pcolor(X, Y, C_mesh);
shading interp;
colormap(hot);
caxis([c_min, c_max]);
axis equal tight;
hold on;
plot(R*cos(theta), R*sin(theta), 'w-', 'LineWidth', 2);
hold off;
cb2 = colorbar;
cb2.Label.String = 'Concentration';
cb2.Label.FontSize = 12;
xlabel('x (cm)', 'FontSize', 13);
ylabel('y (cm)', 'FontSize', 13);
h_title = title(sprintf('%s  |  \\rho = %.2f  |  t = 0.000 s', method_name, rho), 'FontSize', 14);

% Animation
frames_to_show = 200;
step = max(1, floor(size(sol, 2) / frames_to_show));

for j = 1:step:size(sol, 2)
     if ~isvalid(h_pcolor)      % window was closed - stop cleanly
        break
     end
     
    t_now = (j-1) * dt;
    c_now = sol(:, j)';
    C_mesh = repmat(c_now, N_theta, 1);
    
    set(h_pcolor, 'CData', C_mesh);
    set(h_title, 'String', sprintf('%s  |  \\rho = %.2f  |  t = %.3f s', method_name, rho, t_now));
    
    drawnow;
    pause(0.02);
end

fprintf('\n=== Simulation Complete ===\n');
fprintf('Method: %s\n', method_name);
fprintf('D = %.2f, k = %.2f, rho = %.2f\n', D, k, rho);
fprintf('Cell radius: %.2f cm\n', R);
fprintf('Final time: %.3f s\n', t_end);
fprintf('Grid: %d radial points, %d time steps\n', Nr+1, Nt);
