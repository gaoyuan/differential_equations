clear all;close all;clc
%% initial configuration

L = 20; % the length of interval
m = 64; % the number of spatial points 
sd = L / m; % the spatial stepping
v = 0.001; % parameters in the equation

%% The matrices of Laplace, dx, dy operators

% computing Laplace operator LA
n = m * m;
e0 = zeros(n, 1);
e1 = ones(n, 1);
e2 = e1; e4 = e0;

for j = 1:m
    e2(m * j) = 0;
    e4(m * j) = 1;
end

e3(2:n, 1) = e2(1:n-1, 1);
e3(1, 1) = e2(n, 1);
e5(2:n, 1) = e4(1:n-1, 1);
e5(1, 1) = e4(n, 1);

LA = spdiags([e1 e1 e5 e2 -4*e1 e3 e4 e1 e1], [-(n-m) -m -m+1 -1 0 1 m-1 m (n-m)],n,n);
LA = LA / (sd^2);

% computing the dx operator DX
DX = spdiags([e1 -e1  e1 -e1],[-(n-m) -m  m (n-m)],n,n);
DX = DX / sd / 2;

% computing the dy operator DY
DY = spdiags([e5 -e2  e3 -e4],[-m+1 -1  1 m-1],n,n);
DY = DY / sd / 2;

%% solving the system using FFT and ode45

% forming the initial guess
x = linspace(-L/2, L/2, m + 1);
x = x(1:m);
y = linspace(-L/2, L/2, m + 1);
y = y(1:m);

[X,Y] = meshgrid(x,y);
W = exp(-X.^2 - Y.^2/20) + exp(-(X+3*ones(m)).^2 - Y.^2/20);

% generating random initial guess
% W = zeros(m,m);
% for i = 1:10
%     sgn = rand(1);
%     if sgn > 0.5
%         sgn = 1;
%     else
%         sgn = -1;
%     end
%     W = W + sgn * exp(-(X+40*rand(1)*ones(m)).^2/(rand(1)*60) - (Y+40*rand(1)*ones(m)).^2/(rand(1)*60));
% end

W = reshape(W, m*m, 1);

% the K to use in getting psi
k = (2*pi/L)*[0:(m/2-1) (-m/2):-1];
k(1) = 1e-6;

[kX,kY] = meshgrid(k,k);
K = kX.^2 + kY.^2;

% step forward in time using ode45
tspan = 0:0.5:60;
[t,yy] = ode45('rhs', tspan, W, [], K, LA, DX, DY, v, m);

%% Visualizing the result

filename = 'movie.gif';

figure(1); clf;
%set(gcf,'Color','w');
%set(gca,'nextplot','replacechildren','visible','off')
%set(gcf,'Position',[0 0 150 150]);
%set(gca,'Zlim',[-1 1])
solution = reshape(yy(1,:),m,m);
imagesc(x,y,solution)
f = getframe(figure(1));
[im,map] = rgb2ind(f.cdata,256,'nodither');

for i = 2:121
    solution = reshape(yy(i,:),m,m);
    imagesc(x,y,solution)
    %axis('Position',[.1 .1 .75 .8],'Visible','off');
    %subplot(1,2,1)
    %pcolor(X,Y,solution);
    %shading flat
    %axis([-L/2 L/2 -L/2 L/2])
    %axis('Position',[.1 .1 .75 .8],'Visible','off');
    %axis([-L/2 L/2 -L/2 L/2 -1 1 1 100])
    %subplot(1,2,2)
    %surf(X,Y,solution);
    %shading flat
    %colormap(1-0.7*hot)
    %axis([-L/2 L/2 -L/2 L/2 -1 1])
    %set(gcf,'Color','w');
    %set(gca,'Zlim',[-1 1])
    %set(gcf,'Position',[0 0 150 150]);
    f = getframe(figure(1));
    im(:,:,1,i-1) = rgb2ind(f.cdata,map,'nodither');
end

imwrite(uint8(im),map,filename,'LoopCount',inf,'DelayTime',0.05);





