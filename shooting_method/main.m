% Calculate eigenvalues and eigenfunctions using shooting scheme.
% Adjust epsilon due to observation of the boundary value at right.
close all;clear all; clc
K = 1; % the parameter in the equation
tol = 1e-4;
xp = -4:0.1:4; % xspan
epsilon_start = 0.5; % the initial guess of the first(smallest) eigenvalue
A1 = []; A2 = []; % A1: eigenfunctions, A2: eigenvalues
for mode = 1:5 % calculate first five mode
    % The on_off is a switch. In mode 1, the value at right increases
    % as epsilon increases. In mode 2, it's the opposite.
    on_off = (-1)^(mode + 1);
    epsilon = epsilon_start;
    delta = 1; % stepping of epsilon, shrink when overshoot
    for j = 1:1000
        x0 = [1 sqrt(K*xp(1)^2-epsilon)]; % value and launch angle at left
        [t,y] = ode45('shoot',xp,x0,[],epsilon, K); % step forward in time
        
        % if the boundary condition is satisfied at right
        if abs(y(end,2)+y(end,1)*sqrt(K*xp(end)^2-epsilon)) < tol
            norm = trapz(t, y(:,1).*y(:,1)); % trapzoid approx. of area
            A1 = [A1 abs(y(:,1)/sqrt(norm))];
            A2 = [A2 epsilon];
            break
        end
        
        % if undershoot
        if  on_off * (y(end,2)+y(end,1)*sqrt(K*xp(end)^2-epsilon)) > 0
            epsilon = epsilon + delta;
        % if overshoot
        else
            delta = delta / 2; % shrink the step
            epsilon = epsilon - delta;
        end
    end
    epsilon_start = epsilon + 0.5; % the initial guess of next eigenvalue
    plot(t,y(:,1)/sqrt(norm));hold on % plot the eigenfunction
end
save A1.dat A1 -ascii
save A2.dat A2 -ascii