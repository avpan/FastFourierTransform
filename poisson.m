clear all;clc;clf;

Lx = 5;
Ly = 5;
N = 64;
L = 10;
R = 1;
Q = 1.6E-19;
e = 8.85e-12;
A = e/Q;
%boundary conditions
x = linspace(-Lx,Lx,N+1);
x = x(1:N);
%x = removerows(x',N/2+1)';
y = linspace(-Ly,Ly,N+1);
y = y(1:N);
%y = removerows(y',N/2+1)';
% r2 = x.^2 + y.^2;

[X Y] = meshgrid(x,y);

%uniform circular charge where q/e0 = 1;
rho_0 = 5;
rho = zeros(N,N);

for i = 1:N
    for j = 1:N
        r = sqrt(X(i,j)^2 + Y(i,j)^2);
        
        if(X(i,j) == 0 && Y(i,j) == 0)
            rho(i,j) = 1;
        elseif(r < R)
            rho(i,j) = exp((X(i,j)^2 + Y(i,j)^2));
            %rho(i,j) = 1/(pi*(X(i,j)^2 + Y(i,j)^2));
        else
            rho(i,j) = 0;
        end
    end
end

%determine Kn matrices
k = (2*pi/L)*[1-N/2:N/2];
khold = k(N/2+1:N);
k(N/2+1:N) = k(1:N/2);
k(1:N/2) = khold;
k = k(1:N);

[KX KY] = meshgrid(k,k);
K2 = KX.^2 + KY.^2;
K = sqrt(K2);
K2(N,N) = 1;

%find potential
Vn = fft(rho);
V = real(ifft(Vn./K2));

[Ex, Ey] = gradient(V);
% EX = real(ifft(Vn.*KX.*i));
% EY = real(ifft(Vn.*KY.*i));

surf(X,Y,V); colormap('white');
hold on; quiver(X,Y,-Ex,-Ey,10,'black','LineWidth',1.5);
xlabel('X');ylabel('Y');zlabel('Potential (V)');
axis([-Lx,Lx,-Ly,Ly,0,.02]);

