%%%___Forward Difference method for 1D Spherical Diffusion-Reaction___%%%
%
% Solves: dc/dt = D * [d2c/dr2 + (2/r)*dc/dr] + k*c
%
% This is your forward_difference_textbook.m modified for spherical coords.
% Key change: the tridiagonal matrix is NO LONGER constant -- each row
% has coefficients that depend on r(i) because of the (2/r)(dc/dr) term.
%
% input:  spatial interval [rl, rr], time interval [tb, te],
%         number of space steps M, number of time steps N,
%         D = diffusion coeff, IC = initial condition function,
%         BC_R_func = boundary condition at r=R, k = reaction rate
% output: solution w

function w = forward_difference_spherical_nsc(rl, rr, tb, te, M, N, D, IC, BC_R_func, k)

dr = (rr - rl) / M;       % spatial step
dt = (te - tb) / N;       % time step

sigma = D * dt / (dr^2);  % base stability parameter
k_dt = k * dt;             % reaction parameter
fprintf('forward sigma = %f\n', sigma);

m = M - 1;   % number of interior points (excludes r=0 and r=R)
n = N;        % number of time steps

%% ==================== BUILD r-DEPENDENT MATRIX ====================
% Interior grid: r(1)=dr, r(2)=2*dr, ..., r(m)=m*dr
% (We skip r=0 because it's handled separately as the symmetry BC)
r_interior = (1:m)' * dr;

% For each interior point i, the finite difference of the spherical Laplacian:
%
%   D*[c(i+1) - 2c(i) + c(i-1)]/dr^2 + D*(2/r(i))*[c(i+1) - c(i-1)]/(2*dr) + k*c(i)
%
% Collecting coefficients of c(i-1), c(i), c(i+1):
%
%   c(i-1) coefficient: sigma * (1 - 1/i)      where i = r_index = r(i)/dr
%   c(i)   coefficient: 1 - 2*sigma + k_dt
%   c(i+1) coefficient: sigma * (1 + 1/i)
%
% Note: i here is the grid index, and r(i) = i*dr, so 2/(r(i)) * dr/2 = 1/i

% Index vector (i = 1, 2, ..., m corresponding to r = dr, 2dr, ..., m*dr)
idx = (1:m)';

% Coefficients for each interior point
coeff_lower = sigma * (1 - 1./idx);    % c(i-1) coefficient
coeff_diag  = (1 - 2*sigma + k_dt);     % c(i) coefficient (same for all)
coeff_upper = sigma * (1 + 1./idx);    % c(i+1) coefficient

% Build the tridiagonal matrix
a = diag(coeff_diag * ones(m,1));                     % main diagonal
a = a + diag(coeff_upper(1:m-1), 1);                  % upper diagonal
a = a + diag(coeff_lower(2:m), -1);                   % lower diagonal

%% ==================== BOUNDARY CONDITIONS ====================
% r = R (membrane): user-defined function of time
t_vec = tb + (0:n) * dt;
rside = BC_R_func(t_vec);
rside = rside(:)';

% r = 0 (center): symmetry condition dc/dr = 0
% Using L'Hopital: as r->0, (2/r)(dc/dr) -> 2*d2c/dr2
% So Laplacian -> 3*d2c/dr2 at r=0
% This means c(r=0) is updated using:
%   c_new(0) = c(0) + dt * [3*D*(c(1)-c(0))/dr^2 + k*c(0)]
% We handle this INSIDE the time loop, not in the matrix.

%% ==================== INITIAL CONDITIONS ====================
% Interior points
w(:,1) = IC(r_interior)';

% Also track the center point separately
c_center = IC(0);  % concentration at r = 0

%% ==================== TIME STEPPING ====================
for j = 1:n
    % Interior update: w(:,j+1) = a * w(:,j) + boundary contributions
    
    % Boundary contribution vector
    sides = zeros(m, 1);
    
    % First interior point (i=1) gets contribution from r=0 (the center)
    sides(1) = coeff_lower(1) * c_center;
    
    % Last interior point (i=m) gets contribution from r=R
    sides(m) = coeff_upper(m) * rside(j);
    
    % Matrix multiply + boundary terms
    w(:,j+1) = a * w(:,j) + sides;
    
    % Update center point using L'Hopital symmetry condition
    % Laplacian at r=0: 3*D*(c(dr) - c(0))/dr^2
    c_center = c_center + dt * (3 * D * (w(1,j) - c_center) / dr^2 + k * c_center);
end

%% ==================== ASSEMBLE FULL SOLUTION ====================
% Prepend center point, append membrane BC
c_center_history = zeros(1, n+1);
c_center_history(1) = IC(0);
% Recompute center history (we only tracked the final value above)
% For proper output, let's reconstruct it
c_ctr = IC(0);
c_center_history(1) = c_ctr;
for j = 1:n
    c_center_history(j+1) = c_ctr + dt * (3 * D * (w(1,j) - c_ctr) / dr^2 + k * c_ctr);
    c_ctr = c_center_history(j+1);
end

w = [c_center_history; w; rside];

% Coordinate vectors
r_full = (0:M) * dr;
t = (0:N) * dt;

%% ==================== PLOT ====================
surf(r_full, t, w');
xlabel('r (radius)'); ylabel('t (time)');
zlabel('Concentration c(r,t)');
axis([rl rr tb te -1 2]);
view(60, 30);
title({sprintf('spherical reaction-diffusion, Forward Difference No Space Constant'), ...
       sprintf('D=%.2g,  k=%.2g,  \\sigma=%.3g', D, k, sigma)});
cb=colorbar; cb.Label.String='c(r,t)';
view(50, 30);  grid on;
zlim([min(0,min(w(:))) max(w(:))*1.05]);   % data-set, not fixed 
end
