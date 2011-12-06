function rhs = rhs(tspan, w, dummy, K, LA, DX, DY, v, n)
    w = reshape(w, n, n);
    fw = fft2(w); % the transformed omega
    psi = ifft2(- fw ./ K); % the psi at a particular time
    w = reshape(w, n*n, 1);
    psi = reshape(psi, n*n, 1);
    rhs = v*LA*w - (DX*psi).*(DY*w) + (DY*psi).*(DX*w);