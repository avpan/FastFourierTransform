%
% 2D Poisson eqn -(u_xx+u_yy) = f on [0,1]x[0,1], hx = 1/M & hy = 1/N 
%
% this code ASSUMES u=0 on the boundary (Gamma)
%
M = 1024;  hx = 1/M;  x = (0:hx:1)';  lamx = 2*(1-cos(x(2:M)*pi))/(hx^2);
N = 768;  hy = 1/N;  y = (0:hy:1)';  lamy = 2*(1-cos(y(2:N)*pi))/(hy^2);

% Fouriter multiple matrix:  lambda_x(i)+lambda_y(j)

lamx_p_lamy = repmat(lamx,1,N-1) + repmat(lamy',M-1,1);

% exact solution TU with homogeneous boundary value

u = zeros(M+1,N+1);  tu = zeros(M+1,N+1);  f = zeros(M+1,N+1);

for i = 1:M+1
    for j = 1:N+1
          tu(i,j) = sin(3*pi*x(i))*sin(8*pi*y(j));
          f(i,j)  = (3*3+8*8)*pi*pi*tu(i,j);
    end
end

% tell them "what's happening"

format long e
disp('Exact soltuion: u(x,y) = sin(3*pi*x)*sin(8*pi*y)');
disp(['hx = 1/M, M=',num2str(M),'  hy=1/N, N=',num2str(N)]);

% setup  RHS and transform

tic

f = f(2:M,2:N);
fhat =  dst(dst(f)')';

% set u_hat = f_hat divide by Fourier multiples

uhat = fhat ./ (lamx_p_lamy);

% transform back, place in u, then check and print error

u(2:M,2:N) =  dst(dst(uhat)')';
toc 

disp(['Infinty-norm error : ',num2str(max(max(abs(u-tu)))) ]);  


