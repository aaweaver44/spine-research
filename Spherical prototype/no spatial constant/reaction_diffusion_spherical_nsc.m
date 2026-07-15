clear all; close all;

%%%__Initialize parameters, grid, & ICs for SPHERICAL diffusion__%%%
%
% Solves: dc/dt = D/r^2 * d/dr(r^2 * dc/dr) + k*c
%       = D * [d2c/dr2 + (2/r)*dc/dr] + k*c
%
% Domain: r = 0 (center of cell) to r = R (cell membrane)

D = 1;         % diffusion coefficient (cm^2/s)
k = 1;         % reaction rate (1/s), positive = growth
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

%%%__Numerical Solutions__%%%

% Crank-Nicolson (implicit, spherical)
figure;
crank_matrix = crank_nicholson_spherical_nsc(r_start, r_end, t_start, t_end, ...
                                          Nr, Nt, D, IC, BC_R, k);

% Forward Euler (explicit, spherical)
figure;
forward_matrix = forward_difference_spherical_nsc(r_start, r_end, t_start, t_end, ...
                                               Nr, Nt, D, IC, BC_R, k);

%%%__Compare the two methods at final time (line plot)__%%%
figure;
r_full = linspace(r_start, r_end, Nr+1);
plot(r_full, crank_matrix(:,end), 'b-o', 'LineWidth', 1.5, 'DisplayName', 'Crank-Nicolson');
hold on;
plot(r_full, forward_matrix(:,end), 'r--x', 'LineWidth', 1.5, 'DisplayName', 'Forward Difference');
xlabel('r (radius)');
ylabel('Concentration c(r, t_{end})');
title('Comparison at Final Time');
legend('Location', 'best');
grid on;
hold off;

%%%==========================================================%%%
%%%      CIRCULAR CELL VISUALIZATION (Animated)               %%%
%%%==========================================================%%%
%
% Because c only depends on r (spherical symmetry), the concentration
% is the same in every direction. We map c(r) onto a 2D disk by
% creating a polar mesh and assigning c based on each point's distance
% from the center.
%
% Think of it as slicing the sphere in half and looking at the
% cross-section -- every ring at distance r has the same color.

% Choose which solution to animate (change to forward_matrix if desired)
sol = crank_matrix;   
method_name = 'Crank-Nicolson';

% Build polar mesh for the circular cross-section
N_theta = 200;                              % angular resolution
theta = linspace(0, 2*pi, N_theta);         % full circle
[R_mesh, Theta_mesh] = meshgrid(r_full, theta);   % polar grid
X = R_mesh .* cos(Theta_mesh);              % convert to Cartesian x
Y = R_mesh .* sin(Theta_mesh);              % convert to Cartesian y

% Consistent color limits across all plots
c_min = min(sol(:));
c_max = max(sol(:));
if c_min == c_max
    c_max = c_min + 1;  % avoid error if solution is flat
end

% ---- Figure: Snapshot grid (6 time slices) ----
figure('Color', 'w', 'Position', [50, 50, 1200, 700]);
sgtitle(['Cell Cross-Section Over Time (', method_name, ')'], ...
        'FontSize', 16, 'FontWeight', 'bold');

num_snapshots = 6;
snapshot_indices = round(linspace(1, size(sol, 2), num_snapshots));

for s = 1:num_snapshots
    j = snapshot_indices(s);
    t_now = (j-1) * dt;
    
    % Get c(r) at this time step and map onto 2D disk
    c_snapshot = sol(:, j)';                         % row vector of c(r)
    C_mesh = repmat(c_snapshot, N_theta, 1);         % same c at every angle
    
    subplot(2, 3, s);
    pcolor(X, Y, C_mesh);
    shading interp;
    colormap(hot);
    caxis([c_min, c_max]);
    axis equal tight off;
    
    % Draw membrane circle
    hold on;
    plot(R*cos(theta), R*sin(theta), 'w-', 'LineWidth', 1.5);
    hold off;
    
    title(sprintf('t = %.3f s', t_now), 'FontSize', 12);
end

% Shared colorbar
cb = colorbar('Position', [0.93, 0.15, 0.02, 0.7]);
cb.Label.String = 'Concentration';
cb.Label.FontSize = 12;

% ---- Figure: Full animation ----
fig_anim = figure('Color', 'w', 'Position', [100, 100, 800, 700]);

% Initial frame
c_init = sol(:, 1)';
C_mesh = repmat(c_init, N_theta, 1);

h_pcolor = pcolor(X, Y, C_mesh);
shading interp;
colormap(hot);
caxis([c_min, c_max]);
axis equal tight;
hold on;
h_membrane = plot(R*cos(theta), R*sin(theta), 'w-', 'LineWidth', 2);
hold off;
cb2 = colorbar;
cb2.Label.String = 'Concentration';
cb2.Label.FontSize = 12;
xlabel('x (cm)', 'FontSize', 13);
ylabel('y (cm)', 'FontSize', 13);
h_title = title(sprintf('%s  |  t = 0.000 s', method_name), 'FontSize', 14);

% Animate through time
frames_to_show = 200;                                    % total animation frames
step = max(1, floor(size(sol, 2) / frames_to_show));     % skip frames for speed

for j = 1:step:size(sol, 2)
    t_now = (j-1) * dt;
    
    c_now = sol(:, j)';
    C_mesh = repmat(c_now, N_theta, 1);
    
    set(h_pcolor, 'CData', C_mesh);
    set(h_title, 'String', sprintf('%s  |  t = %.3f s', method_name, t_now));
    
    drawnow;
    pause(0.02);
end

fprintf('\n=== Simulation Complete ===\n');
fprintf('Method: %s\n', method_name);
fprintf('Cell radius: %.2f cm\n', R);
fprintf('Final time: %.3f s\n', t_end);
fprintf('Grid: %d radial points, %d time steps\n', Nr+1, Nt);
