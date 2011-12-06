function rhs = rhs(tspan, rhs0, dummy, K, c, n)
    
    rhs0 = reshape(rhs0, n, n, n);
    irhs0 = ifftn(rhs0);
    rhs = -1i*(0.5*K.*rhs0 + fftn(abs(irhs0).^2.*irhs0 - c.*irhs0));
    rhs = reshape(rhs, n^3, 1);
    