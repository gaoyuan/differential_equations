function rhs = shoot(xp,x,dummy,epsilon, K)
% x(1) is the solution phi.
% Let x(1)' = x(2).
% Then the equation tells us:
% x(2)' = x(1)'' = (Kx^2 - epsilon)x(1).
    rhs = [x(2);x(1)*(K * xp^2-epsilon)];