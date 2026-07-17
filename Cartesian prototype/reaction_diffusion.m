clear all; close all;
    %%%__Initialize parameters, grid, & ICs__%%%
D = 1; % diffusion coefficient (cm^2/s)
k_n = 1; % nutrients feeding the cell (1/s) where negative = consumption, positive = influx
L = 1; % domain length (cm)

x_start = 0; % start pos
x_end = x_start+L; % end pos

dx = .1; % spatial step size (cm)
Nx = round(L/dx); % number of spatial steps
x = linspace(x_start, x_end, Nx+1); % matrix of spatial grid

t_start = 0;
t_end = 1; % sec; end simulation after this many seconds have elapsed
dt_req = 0.004; % requested time step size (sec) 
Nt = round(t_end/dt_req); % number of time steps 
dt = t_end/Nt;  % actual time step size
t_matrix = linspace(t_start, t_end, Nt+1); % matrix of time steps
fprintf('dt: requested %.6f -> effective %.6f  (Nt = %d),  sigma = %.4f\n', ...
        dt_req, dt, Nt, D*dt/dx^2);

c_matrix = zeros(Nx, Nt); % Initialize concentration matrix
% Initial Conditions
c0 = 1; % initial concentration
lbc = 0; % left boundary condition
rbc = 0; % right boundary condition
IC = @(x)sin(2*pi*x).^2; % initial condition global fcn ; u(x,0)=c0 (same size as inpout)
BC_L = @(t)lbc*t; % left BC global fcn ; u(0,t)=lbc
BC_R = @(t)rbc*t; % right BC global fcn ; u(1,t)=rbc

%%%__Numerical solutions__%%%

% Backward Euler (implicit)

% Crank-Nicolson / trapezoid
figure;
crank_matrix = crank_nicholson_textbook(x_start, x_end, t_start, t_end, Nx, Nt, D, IC, BC_L, BC_R, k_n); % contains a matrix of u(x,t) (spatial by time)

% Runge-Kutta

% Forward Euler (explicit)
figure;
forward_matrix = forward_difference_textbook(x_start, x_end, t_start, t_end, Nx, Nt, D, IC, BC_L, BC_R ,k_n); % contains a matrix of u(x,t) (spatial by time)


