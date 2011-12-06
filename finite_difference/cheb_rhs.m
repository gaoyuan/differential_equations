function rhs = cheb_rhs(t, s, dummy, L, N, D1, D2, beta)
    u = s(1:N*N);
    v = s(N*N+1:2*N*N);
    A = u.^2 + v.^2;
    ut = (ones(N*N,1) - A).*u + beta * A .* v + D1 * L * u / 100;
    vt = - beta * A .* u + (ones(N*N,1) - A).*v  + D2 * L * v / 100;
    rhs = [ut;vt];