function rhs = shoot_nonlinear(xp,x,dummy,gama,epsilon, K)
    rhs = [x(2);gama*x(1).^3 + x(1)*(K*xp^2-epsilon)];