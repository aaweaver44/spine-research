%%%___Crank-Nicolson method for 1D Spherical Diffusion-Reaction___%%%
%
% Solves: dc/dt = D * [d2c/dr2 + (2/r)*dc/dr] + k*c
%
% This is your crank_nicholson_textbook.m modified for spherical coords.
% Key change: the tridiagonal matrices A (implicit) and B (explicit) now
% have r-dependent coefficients from the (2/r)(dc/dr) spherical term.
%
% input:  spatial interval [rl, rr], time interval [tb, te],
%         number of space steps M, number of time steps N,
%         D = diffusion coeff, IC = initial condition function,
%         BC_R_func = boundary condition at r=R, k = reaction rate
% output: solution w

function w = crank_nicholson_spherical(rl, rr, tb, te, M, N, D, IC, BC_R_func, k)

dr = (rr - rl) / M;       % spatial step
dt = (te - tb) / N;       % time step

sigma = D * dt / (dr^2);  % base stability parameter
rho = k * dt;             % reaction parameter
fprintf('crank sigma = %f\n', sigma);

m = M - 1;   % number of interior points
n = N;        % number of time steps

%% ==================== BUILD r-DEPENDENT MATRICES ====================
% Interior grid: r(i) = i*dr for i = 1, 2, ..., m
% Index vector
idx = (1:m)';

% Spherical correction factor at each grid point: 1/i comes from
%   (2/r_i) * (dr / 2) = (2/(i*dr)) * (dr/2) = 1/i
%
% The stencil coefficients for the spherical Laplacian:
%   c(i-1): sigma * (1 - 1/i)
%   c(i):   -2*sigma
%   c(i+1): sigma * (1 + 1/i)
%
% Crank-Nicolson averages the spatial operator between time levels j and j+1:
%
% IMPLICIT side (time j+1):  [2I - dt*L] * w^(j+1) = ...
%   c(i-1): -sigma * (1 - 1/i)
%   c(i):    2 + 2*sigma - rho
%   c(i+1): -sigma * (1 + 1/i)
%
% EXPLICIT side (time j):    ... = [2I + dt*L] * w^(j) + BC terms
%   c(i-1):  sigma * (1 - 1/i)
%   c(i):    2 - 2*sigma + rho
%   c(i+1):  sigma * (1 + 1/i)

% Coefficients that vary with r
lower_coeff = sigma * (1 - 1./idx);   % sub-diagonal
upper_coeff = sigma * (1 + 1./idx);   % super-diagonal

%% ---- Implicit matrix A (left-hand side) ----
A = diag((2 + 2*sigma - rho) * ones(m,1));           % main diagonal
A = A + diag(-upper_coeff(1:m-1), 1);                 % upper diagonal
A = A + diag(-lower_coeff(2:m), -1);                  % lower diagonal

%% ---- Explicit matrix B (right-hand side) ----
B = diag((2 - 2*sigma + rho) * ones(m,1));            % main diagonal
B = B + diag(upper_coeff(1:m-1), 1);                  % upper diagonal
B = B + diag(lower_coeff(2:m), -1);                   % lower diagonal

%% ==================== BOUNDARY CONDITIONS ====================
% r = R (membrane): user-defined function of time
t_vec = tb + (0:n) * dt;
rside = BC_R_func(t_vec);
rside = rside(:)';

% r = 0 (center): symmetry condition dc/dr = 0
% Using L'Hopital: Laplacian at r=0 -> 3*d2c/dr2
% We handle the center point separately in the time loop.

%% ==================== INITIAL CONDITIONS ====================
r_interior = (1:m)' * dr;
w(:,1) = IC(r_interior)';

% Track center point separately
c_center = IC(0);
c_center_history = zeros(1, n+1);
c_center_history(1) = c_center;

%% ==================== TIME STEPPING ====================
for j = 1:n
    % Boundary contribution vector for Crank-Nicolson
    % Both time levels j and j+1 contribute
    sides = zeros(m, 1);
    
    % First interior point (i=1): contribution from center (r=0)
    % The center value at time j is known; at time j+1 we use current estimate
    % For simplicity, we use a semi-implicit approach for the center coupling
    sides(1) = lower_coeff(1) * c_center;  % explicit side contribution
    
    % Last interior point (i=m): contribution from membrane (r=R)
    % Crank-Nicolson uses average of both time levels
    sides(m) = upper_coeff(m) * (rside(j) + rside(j+1));
    
    % Also add the center contribution to the implicit side correction
    sides_implicit_correction = zeros(m, 1);
    sides_implicit_correction(1) = lower_coeff(1) * c_center;  % will update after
    
    % Right-hand side
    rhs = B * w(:,j) + sides + sides_implicit_correction;
    
    % Solve the implicit system
    w(:,j+1) = A \ rhs;
    
    % Update center point using L'Hopital symmetry (Crank-Nicolson style)
    % At r=0: dc/dt = 3*D*(d2c/dr2) + k*c
    % CN average: c_new = c_old + dt/2 * [L(c_old) + L(c_new)]
    % where L(c) at center = 3*D*(c(dr) - c(0))/dr^2 + k*c(0)
    % Using the new w(1,j+1) for the implicit part:
    L_old = 3 * D * (w(1,j) - c_center) / dr^2 + k * c_center;
    
    % Predict new center (need to solve implicitly)
    % c_new = c_old + dt/2 * [L_old + 3*D*(w(1,j+1) - c_new)/dr^2 + k*c_new]
    % c_new * [1 + dt/2*(3*D/dr^2 - k)] = c_old + dt/2*L_old + dt/2*3*D*w(1,j+1)/dr^2
    coeff_implicit = 1 + dt/2 * (3*D/dr^2 - k);
    rhs_center = c_center + dt/2 * L_old + dt/2 * 3*D * w(1,j+1) / dr^2;
    c_center = rhs_center / coeff_implicit;
    c_center_history(j+1) = c_center;
end

%% ==================== ASSEMBLE FULL SOLUTION ====================
w = [c_center_history; w; rside];

% Coordinate vectors
r_full = (0:M) * dr;
t = (0:N) * dt;

%% ==================== PLOT ====================
surf(r_full, t, w');
xlabel('r (radius)'); ylabel('t (time)');
zlabel('Concentration c(r,t)');
axis([rl rr tb te -1 2]);
title('Crank-Nicolson - Spherical Coordinates');
colorbar;

end
