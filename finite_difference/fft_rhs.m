function rhs = fft_rhs(t, s, dummy, N, K, D1, D2, beta)
    U = reshape(s(1:N*N), N, N);
    V = reshape(s(N*N+1:2*N*N), N, N);
    U = ifft2(U); V = ifft2(V);
    A = U.^2 + V.^2;
    lambda = ones(N) - A;
    omega = - beta * A;
    Ut = fft2(lambda .* U - omega .* V) - D1 * K .* fft2(U);
    Vt = fft2(omega .* U + lambda .* V) - D2 * K .* fft2(V);
    Ut = reshape(Ut, N*N, 1);
    Vt = reshape(Vt, N*N, 1);
    rhs = [Ut;Vt];