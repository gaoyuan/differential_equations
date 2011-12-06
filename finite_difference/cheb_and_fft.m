clear all; close all; clc

L = 20;
N = 64;

%% periodic boundaries with fft
x = linspace(-L/2, L/2, N + 1);
x = x(1:N);
y = linspace(-L/2, L/2, N + 1);
y = y(1:N);
[X, Y] = meshgrid(x, y);

m = 1; % number of spirals
beta = 1; % coeefficient for omega (physical meaning?)
D1 = 0.1; % diffusion coefficient
D2 = 0.1; % diffusion coefficient
u = tanh(sqrt(X.^2+Y.^2)).*cos(m*angle(X+i*Y)-(sqrt(X.^2+Y.^2)));
v = tanh(sqrt(X.^2+Y.^2)).*sin(m*angle(X+i*Y)-(sqrt(X.^2+Y.^2)));
U = reshape(fft2(u), N*N, 1);
V = reshape(fft2(v), N*N, 1);
sol0 = [U;V]; % initial solution in the Fourier domain

k = (2*pi/L)*[0:(N/2-1) (-N/2):-1];
[kX,kY] = meshgrid(k,k);
K = kX.^2 + kY.^2;

tspan = 0:.5:4;
% solution in the Fourier domain
[t, sol] = ode45('fft_rhs', tspan, sol0, [], N, K, D1, D2, beta);

% plotting the result for periodic boundary solutions
% for i = 1:9
%     U = ifft2(reshape(sol(i,1:N*N), N, N));
%     V = ifft2(reshape(sol(i,N*N+1:2*N*N), N, N));
%     pcolor(X,Y,U); shading interp
%     pause(0.5)
% end

%% no-flux boundaries with cheb
N = 31;
[D, x] = cheb(N - 1);
DD = D*D;
DD(1, :) = zeros(1, N); % no-flux
DD(end, :) = zeros(1, N); % no-flux
% m = 1; % number of spirals
% D1 = 0.1; % diffusion coefficient
% D2 = 0.1; % diffusion coefficient
% beta = 1;
% tspan = 0:.5:4;
[X,Y] = meshgrid(x*L/2,x*L/2);
u = tanh(sqrt(X.^2+Y.^2)).*cos(m*angle(X+i*Y)-(sqrt(X.^2+Y.^2)));
v = tanh(sqrt(X.^2+Y.^2)).*sin(m*angle(X+i*Y)-(sqrt(X.^2+Y.^2)));
sol0 = [reshape(u, N*N, 1);reshape(v, N*N, 1)];
I = eye(length(DD));
L = kron(I,DD) + kron(DD,I); % 2D Laplacian
[t, sol] = ode45('cheb_rhs', tspan, sol0, [], L, N, D1, D2, beta);

% To visualize...
% size(sol)
% U = reshape(sol(9,1:N*N), N, N);
% V = reshape(sol(9,N*N+1:2*N*N), N, N);
% pcolor(X,Y,U); shading interp

