%%%___Forward Difference method for 1D Spherical Diffusion-Reaction___%%%
%
% Solves: dc/dt = D * [d2c/dr2 + (2/r)*dc/dr] + k*c + rho
%
% rho is a constant spatial source/sink term (units: concentration/time)
%
% input:  spatial interval [rl, rr], time interval [tb, te],
%         number of space steps M, number of time steps N,
%         D = diffusion coeff, IC = initial condition function,
%         BC_R_func = boundary condition at r=R, k = reaction rate,
%         rho = constant source term
% output: solution w

function w = forward_difference_spherical(rl, rr, tb, te, M, N, D, IC, BC_R_func, k, rho)

dr = (rr - rl) / M;       % spatial step
dt = (te - tb) / N;       % time step

sigma = D * dt / (dr^2);  % base stability parameter
rho_k = k * dt;           % reaction parameter (renamed to avoid conflict with rho)
rho_dt = rho * dt;         % source term scaled by time step
fprintf('forward sigma = %f\n', sigma);

m = M - 1;   % number of interior points
n = N;        % number of time steps

%% ==================== BUILD r-DEPENDENT MATRIX ====================
idx = (1:m)';

% Tridiagonal coefficients from spherical Laplacian
coeff_lower = sigma * (1 - 1./idx);    % c(i-1) coefficient
coeff_diag  = (1 - 2*sigma + rho_k);   % c(i) coefficient
coeff_upper = sigma * (1 + 1./idx);    % c(i+1) coefficient

% Build the tridiagonal matrix
a = diag(coeff_diag * ones(m,1));
a = a + diag(coeff_upper(1:m-1), 1);
a = a + diag(coeff_lower(2:m), -1);

%% ==================== BOUNDARY CONDITIONS ====================
t_vec = tb + (0:n) * dt;
rside = BC_R_func(t_vec);
rside = rside(:)';

%% ==================== INITIAL CONDITIONS ====================
r_interior = (1:m)' * dr;
w(:,1) = IC(r_interior)';
c_center = IC(0);

%% ==================== TIME STEPPING ====================
for j = 1:n
    % Boundary contribution vector
    sides = zeros(m, 1);
    sides(1) = coeff_lower(1) * c_center;
    sides(m) = coeff_upper(m) * rside(j);
    
    % Constant source term vector (rho * dt added at every interior point)
    source = rho_dt * ones(m, 1);
    
    % Matrix multiply + boundary terms + constant source
    w(:,j+1) = a * w(:,j) + sides + source;
    
    % Update center point using L'Hopital symmetry condition
    % Laplacian at r=0: 3*D*(c(dr) - c(0))/dr^2
    c_center = c_center + dt * (3 * D * (w(1,j) - c_center) / dr^2 ...
               + k * c_center + rho);
end

%% ==================== ASSEMBLE FULL SOLUTION ====================
c_ctr = IC(0);
c_center_history = zeros(1, n+1);
c_center_history(1) = c_ctr;
for j = 1:n
    c_center_history(j+1) = c_ctr + dt * (3 * D * (w(1,j) - c_ctr) / dr^2 ...
                             + k * c_ctr + rho);
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
title(sprintf('Forward Difference - Spherical  (\\rho = %.2f)', rho));
colorbar;

end
