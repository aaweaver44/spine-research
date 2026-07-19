%%%___Crank-Nicolson method for 1D Spherical Diffusion-Reaction___%%%
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

function w = crank_nicholson_spherical_p(rl, rr, tb, te, M, N, D, IC, BC_R_func, k, rho)

dr = (rr - rl) / M;       % spatial step
dt = (te - tb) / N;       % time step

sigma = D * dt / (dr^2);  % base stability parameter
k_dt = k * dt;           % reaction parameter (renamed to avoid conflict with rho)
rho_dt = rho * dt;         % source term scaled by time step
fprintf('crank sigma = %f\n', sigma);

m = M - 1;   % number of interior points
n = N;        % number of time steps

%% ==================== BUILD r-DEPENDENT MATRICES ====================
idx = (1:m)';

% Spherical correction coefficients
lower_coeff = sigma * (1 - 1./idx);
upper_coeff = sigma * (1 + 1./idx);

%% ---- Implicit matrix A (left-hand side) ----
A = diag((2 + 2*sigma - k_dt) * ones(m,1));
A = A + diag(-upper_coeff(1:m-1), 1);
A = A + diag(-lower_coeff(2:m), -1);

%% ---- Explicit matrix B (right-hand side) ----
B = diag((2 - 2*sigma + k_dt) * ones(m,1));
B = B + diag(upper_coeff(1:m-1), 1);
B = B + diag(lower_coeff(2:m), -1);

%% ==================== BOUNDARY CONDITIONS ====================
t_vec = tb + (0:n) * dt;
rside = BC_R_func(t_vec);
rside = rside(:)';

%% ==================== INITIAL CONDITIONS ====================
r_interior = (1:m)' * dr;
w(:,1) = IC(r_interior)';

c_center = IC(0);
c_center_history = zeros(1, n+1);
c_center_history(1) = c_center;

%% ==================== TIME STEPPING ====================
for j = 1:n
    % Boundary contribution vector
    sides = zeros(m, 1);
    sides(1) = lower_coeff(1) * c_center;
    sides(m) = upper_coeff(m) * (rside(j) + rside(j+1));
    
    % Implicit correction for center coupling
    sides_implicit_correction = zeros(m, 1);
    sides_implicit_correction(1) = lower_coeff(1) * c_center;
    
    % Constant source: Crank-Nicolson averages over both time levels,
    % but since rho is constant, the average is just rho.
    % Multiply by 2*dt to match the CN scaling (both sides multiplied by 2)
    source = rho_dt * 2 * ones(m, 1);
    
    % Right-hand side
    rhs = B * w(:,j) + sides + sides_implicit_correction + source;
    
    % Solve the implicit system
    w(:,j+1) = A \ rhs;
    
    % Update center point (Crank-Nicolson style)
    L_old = 3 * D * (w(1,j) - c_center) / dr^2 + k * c_center + rho;
    
    coeff_implicit = 1 + dt/2 * (3*D/dr^2 - k);
    rhs_center = c_center + dt/2 * L_old + dt/2 * (3*D * w(1,j+1) / dr^2 + rho);
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
h = surf(r_full, t, w');
h.FaceColor = 'interp';    % smooth color gradient
h.EdgeColor = 'k';         % but keep the black grid lines
h.EdgeAlpha = 0.3;         % faint grid lines
xlabel('r (radius)'); ylabel('t (time,s)'); zlabel('Concentration c(r,t)');
title({'spherical reaction-diffusion, Crank Nicholson (with source)', ...
       sprintf('D=%.2g,  k=%.2g,  rho=%.2g,  \\sigma=%.3g', D, k, rho, sigma)});
cb=colorbar; cb.Label.String='c(r,t)';
view(50, 30);  grid on;
zlim([min(0,min(w(:))) max(w(:))*1.05]);

end
