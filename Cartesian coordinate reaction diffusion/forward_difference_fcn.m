%%%___Program 8.1 Forward Difference method for heat equation___%%%
% input: space interval [xl,xr], time interval [yb,yt],
%        number of space steps M, number of time steps N
% output: solution w
% Example usage: w=heatfßd(0,1,0,1,10,250)
% 8.1: the step size k = 0.004 corresponds to N = 250 steps in time from t = 0 to t = 1.
function w=forward_difference_textbook(xl,xr,yb,yt,M,N,D,IC,BC_L,BC_R,k)
dx=(xr-xl)/M; % h
dt=(yt-yb)/N; % k

sigma=D*dt/(dx*dx);
rho=k*dt;
fprintf('forward sigma = %f\n', sigma);

m=M-1; 
n=N;  

a=diag((1-2*sigma+rho)*ones(m,1))+diag(sigma*ones(m-1,1),1);
a=a+diag(sigma*ones(m-1,1),-1);                        % define matrix a

t_vec = yb+(0:n)*dt;
x_interior = xl+(1:m)*dx;
lside=BC_L(t_vec); 
lside = lside(:)'; % ensure it is a row
rside=BC_R(t_vec);
rside = rside(:)'; % ensure it is a row
w(:,1)=IC(x_interior)';                                 % initial conditions

for j=1:n
  w(:,j+1)=a*w(:,j)+sigma*[lside(j);
  zeros(m-2,1);rside(j)];
end

w=[lside;w;rside];                                     % attach boundary conds
x=(0:m+1)*dx;
t=(0:n)*dt;

% mesh(x,t,w');
surf(x,t,w');
axis([xl xr yb yt -1 2])                  % 3-D plot of solution w
view(60,30);
title('Forward Difference');
% Use dot notation in functions f, l, and r

% % function u=f(x)
% % u=sin(2*pi*x).^2;
% % 
% % function u=l(t)
% % u=0*t;
% % 
% % function u=r(t)
% % u=0*t;