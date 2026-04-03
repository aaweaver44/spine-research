clear all; close all;
    %%%__Initialize parameters, grid, & ICs__%%%
D = 1; % diffusion coefficient (cm^2/s)
k = 1; % nutrients feeding the cell (1/s) where negative = consumption, positive = influx
L = 1; % domain length (cm)
dx = .1; % spatial step size (cm)
    % good value for crank: 0.1
    % good value for forward: 0.1
dt = 0.003; % time step size (sec) 
    % good value for crank: 0.0055
    % good value for forward: 0.003
t_end = 1; % sec; end simulation after this many seconds have elapsed

x_start = 0; % start pos
x_end = x_start+L; % end pos

Nx = round(L/dx); % number of spatial steps
x = linspace(x_start, x_end, Nx+1); % matrix of spatial grid

t_start = 0;
Nt = round(t_end/dt); % number of time steps 
t_matrix = linspace(t_start, t_end, Nt+1); % matrix of time steps

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
crank_matrix = crank_nicholson_textbook(x_start, x_end, t_start, t_end, Nx, Nt, D, IC, BC_L, BC_R, k); % contains a matrix of u(x,t) (spatial by time)

% Runge-Kutta

% Forward Euler (explicit)
figure;
forward_matrix = forward_difference_textbook(x_start, x_end, t_start, t_end, Nx, Nt, D, IC, BC_L, BC_R ,k); % contains a matrix of u(x,t) (spatial by time)


