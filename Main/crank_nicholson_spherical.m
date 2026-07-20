function w=crank_nicholson_spherical(rl,rr,yb,yt,M,N,D,w0,BC_L,BC_R,k,P,dr,r_vec,inner_bc_select)
% w0 : initial condition vector for this time step (length m)
% rl, rr : spatial domain (radius, um)
% yb, yt : time domain (s)
% M : number of spatial divisions
% N : number of time divisions 
% D : diffusion coefficient (um^2/s)
% w0 : 
% k : reaction coefficient (1/s)
% P : constant source term (uM/s)
% dr : spatial step size (um)
% r_vec : radial grid vector (um)
% inner_bc_select : choose how the function handles the inner boundary condition

% Step size
d=(yt-yb)/N; % k
% stability parameter 
sigma=D*d/(dr*dr); % sigma = D*dt/dr^2
rho = k*d;

% number of interior points
m=M-1; 
n=N;
r_interior = r_vec(2:m+1);   % interior radial points, excluding boundaries

% r-dependent off-diagonal coefficients
sigma_minus = sigma * (1 - dr./r_interior);  % sigma_i^-
sigma_plus  = sigma * (1 + dr./r_interior);  % sigma_i^+

% Handle inner bc calculation per selection
if strcmp(inner_bc_select, 'symmetry')       % r = 0 center: L'Hopital symmetry, no left neighbor
    sigma_minus(1) = 0;
    sigma_plus(1)  = 2*sigma;
elseif strcmp(inner_bc_select, 'dirichlet')  % interface: keep r-dependent coeffs so BC_L applies
else
    error('inner_bc_select must be ''symmetry'' or ''dirichlet''');
end

% Implicit side (j terms, LHS)
a = diag((2+2*sigma-rho)*ones(m,1)) ...
  + diag(-sigma_plus(1:m-1),  1) ...
  + diag(-sigma_minus(2:m),  -1);

% Explicit side (j-1 terms, RHS)
b = diag((2-2*sigma+rho)*ones(m,1)) ...
  + diag( sigma_plus(1:m-1),  1) ...
  + diag( sigma_minus(2:m),  -1);

% Boundary and Initial conditions
t_vec = yb+(0:n)*d;
lside=BC_L(t_vec); 
lside = lside(:)'; % ensure it is a row
rside=BC_R(t_vec);
rside = rside(:)'; % ensure it is a row

w(:,1) = w0;            % initial condition passed in from main loop

%%%___MATH TIME___%%%
for j=1:n
  % Set Aside boundary Conditions
  sides=[sigma_minus(1)*(lside(j)+lside(j+1));
  zeros(m-2,1);
  sigma_plus(end)*(rside(j)+rside(j+1))];
% Solution with 2*P*d term added to RHS at every time step
  w(:,j+1)=a\(b*w(:,j)+sides+2*P*d*ones(m,1));
end

% return filled matrix with BCs attached
w=[lside;w;rside];