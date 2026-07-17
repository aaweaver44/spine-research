%%%___Program 8.2 from textbook___%%%

% Program 8.2 Crank-Nicolson method
% input: space interval [xl,xr], time interval [yb,yt],
%        number of space steps M, number of time steps N
% output: solution w
% Example usage: w=crank(0,1,0,1,10,10)

function w=crank_nicholson_textbook(xl,xr,yb,yt,M,N,D,IC,BC_L,BC_R, k_n)
% in main file, passing it: x_start, x_end, t_start, t_end, Nx, Nt, IC, BC_L, BC_R
% xl, xr : spatial domain [xl,xr]
% yb, yt : time domain [yb,yt]
% M : number of spatial divisions
% N : number of time divisions
% c : diffusion coefficient

% Step sizes 
dx=(xr-xl)/M; % dx
dt=(yt-yb)/N; % dt

% stability parameter - model is best when sigma is around 1
% Also what we are using for the finite difference analysis
sigma=D*dt/(dx*dx); 
rho = k_n*dt;
fprintf('crank sigma = %f\n', sigma);

% number of interior points
m=M-1; 
n=N;

% Implicit side 
a=diag(2+2*sigma*ones(m,1))+diag(-sigma*ones(m-1,1),1);
a=a+diag(-sigma*ones(m-1,1),-1);       % define tridiagonal matrix a

% Explicit side
b=diag(2-2*sigma*ones(m,1))+diag(sigma*ones(m-1,1),1);
b=b+diag(sigma*ones(m-1,1),-1);        % define tridiagonal matrix b

% Boundary and Initial conditions
% % lside=l(yb+(0:n)*k); 
% % rside=r(yb+(0:n)*k);
% % w(:,1)=f(xl+(1:m)*h)';                 % initial conditions
t_vec = yb+(0:n)*dt;
x_interior = xl+(1:m)*dx;
lside=BC_L(t_vec); 
lside = lside(:)'; % ensure it is a row
rside=BC_R(t_vec);
rside = rside(:)'; % ensure it is a row
w(:,1)=IC(x_interior)';  

%%%___MATH TIME___%%%
for j=1:n
  % Set Aside boundary Conditions
  sides=[lside(j)+lside(j+1);
  zeros(m-2,1);
  rside(j)+rside(j+1)];
  % Solution
  w(:,j+1)=a\(b*w(:,j)+sigma*sides);
end

% Fill boundary conditions
w=[lside;w;rside];

% Coordinate vectors for plotting
x=xl+(0:M)*dx;
t=yb+(0:N)*dt;

% 3D surface plot
% mesh(x,t,w');
surf(x, t, w');
shading interp;
xlabel('x  (cm)');
ylabel('t  (s)');
zlabel('c(x,t)');
title({'Crank-Nicolson  |  1D Cartesian diffusion', ...
       sprintf('D = %.3g cm^2/s,  \\Deltax = %.3g cm,  \\Deltat = %.4g s,  \\sigma = %.3g  (unconditionally stable)', ...
               D, dx, dt, sigma), ...
       sprintf('max|c| = %.3g', max(abs(w(:))))});
cb = colorbar;  cb.Label.String = 'c(x,t)';
view(60, 30);  grid on;
axis([xl xr yb yt min(-0.1, min(w(:))) max(1.1, max(w(:)))]);
%axis([xl xr yb yt -1 2])

% BC and IC functions
% % function u=f(x)
% % u=sin(2*pi*x).^2;
% % 
% % function u=l(t)
% % u=0*t;
% % 
% % function u=r(t)
% % u=0*t;