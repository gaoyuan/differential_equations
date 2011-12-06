% calculate the eigenvalue and eigenfunctions using finite difference
clear all;close all;clc
d = 0.1; % mesh size
L = 4; % interval
xspan = -L:d:L; % mesh
N = length(xspan); % number of mesh points
dim = N - 2; % dimension of the solution matrix
A = zeros(dim,dim); % the interior points

% --- boundary conditions ---
A(1,1) = xspan(2) * xspan(2) + 2 / 3 / d^2;
A(1,2) = - 2 / 3 / d^2;
A(dim, dim - 1) = - 2 / 3 / d^2;
A(dim, dim) = (L - d)^2 + 2 / 3 / d^2;

% --- finite difference scheme of interior points ---
for i = 2:dim - 1 % the second line to last but one line
    A(i,i-1) = - 1 / d^2;
    A(i,i) = 2 / d^2 + xspan(i+1) * xspan(i+1);
    A(i,i+1) = - 1 / d^2;
end

% --- calculating eigenvalues and eigenfunctions ---
[psi, epsilon] = eigs(A,5,'sm'); % first five mode
epsilon = diag(epsilon)';
[epsilon, index] = sort(epsilon); % sort eigenvalue from small to large
psi = psi(:,index); % sort the corresponding eigenvalues

% --- adding the boundary points ---
psi = [(4 * psi(1,:) - psi(2,:))/3;psi];
psi = [psi;(4 * psi(end,:) - psi(end - 1,:))/3];

% --- normalize and plot ---
hold on
color = ['k','b','r','y','m'];
for i = 1:5
    norm = trapz(xspan, psi(:,i).*psi(:,i));
    psi(:,i) = abs(psi(:,i))/sqrt(norm);
    plot(xspan, psi(:,i), color(i)); % plot the abs of eigenfunctions
end