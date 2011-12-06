%%%%%%%%%%%%%
% problem c %
%%%%%%%%%%%%%
tol = 1e-4;
K = 1; % parameter in the equation
xp = -2:0.1:2; % xspan
A5 = [];A6 = [];A7 = [];A8 = [];

% in this case, we have to adjust two parameters: A and epsilon
for gama = [-0.05,0.05] % two parameters
    epsilon_start = 0.5;
    for mode = 1:2 % first two modes
        A = 0.5;
        dA = 1;
        on_off = (-1)^(mode + 1); % the switch
        for i = 1:1000
            delta = 1;
            epsilon = epsilon_start;
            for j = 1:1000
                x0 = [A sqrt(xp(1)^2-epsilon)*A];
                [t,y] = ode45('shoot_nonlinear',xp,x0,[],gama,epsilon, K);
                
                % if boundary condition is satisfied
                if abs(y(end,2)+y(end,1)*sqrt(xp(end)^2-epsilon)) < tol
                    norm1 = trapz(t, y(:,1).*y(:,1));
                    break
                end
                
                % if undershoot
                if on_off * (y(end,2)+y(end,1)*sqrt(xp(end)^2-epsilon)) > 0
                    epsilon = epsilon + delta;
                else
                    delta = delta / 2;
                    epsilon = epsilon - delta;
                end
            end
            
            % if the norm is correct
            if abs(norm1 - 1) < tol
                if gama == 0.05
                    A5 = [A5 abs(y(:,1))/sqrt(norm1)];
                    A6 = [A6 epsilon];
                else
                    A7 = [A7 abs(y(:,1))/sqrt(norm1)];
                    A8 = [A8 epsilon];
                end
                break
            end
            
            % if norm < 1, amplify the launch angle
            if norm1 < 1
                A = A + dA;
            else
                dA = dA / 2;
                A = A - dA;
            end
        end
        epsilon_start = epsilon_start + 1; % initial guess for next egv.
    end
end