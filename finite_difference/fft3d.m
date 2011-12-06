clear all; close all; clc

% Gross-Pitaevskii system (nonlinear Schrodinger with potential) in 3D

A1 = -1; A2 = -1; A3 = -1; B1 = 1; B2 = 1; B3 = 1; % parameters in the equation
n = 16; % the number of points in each direction (X,Y,Z)
L = 2*pi; % the interval length of solution in each direction

% periodic boundaries, so (n + 1) ...
x = linspace(-L/2, L/2, n + 1);
y = linspace(-L/2, L/2, n + 1);
z = linspace(-L/2, L/2, n + 1);
x = x(1:n); y = y(1:n); z = z(1:n);

[X, Y, Z] = meshgrid(x, y, z);
psi_0 = cos(X).*cos(Y).*cos(Z); % the initial solution
psi_1 = sin(X).*sin(Y).*sin(Z); % the initial solution 2

% construct the constant that is invarient during advencing in time
iden = ones(n,n,n);
c = (A1 * sin(X).^2 + iden * B1).*(A2 * sin(Y).^2 + iden * B2).*(A3 * sin(Z).^2 + iden * B3);

% constructing kx^2 + ky^2 + kz^2
k = (2*pi/L)*[0:(n/2-1) (-n/2):-1]; % the shifted space
[kX, kY, kZ] = meshgrid(k,k,k);
K = kX.^2 + kY.^2 + kZ.^2; 

tspan = 0:0.5:4; % stepping in time
[t,sol] = ode45('fft3d_rhs', tspan, fftn(psi_0), [], K, c, n);
[t,sol1] = ode45('fft3d_rhs', tspan, fftn(psi_1), [], K, c, n);