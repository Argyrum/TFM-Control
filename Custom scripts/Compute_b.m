function [b, isStable, Kn, Jn] = Compute_b(A, B, C, c, wn, verbose)
% COMPUTE_B [b, isStable, Kn, Jn] - Compute relay hysteresys for equilibrium
%   Inverse of FindWnPrecise

    arguments
        A
        B
        C
        c
        wn
        verbose = true
    end

    Jn = lprsmatr(A, B, C, wn);
    b = -imag(Jn)*4*c/pi;

    if TestOrbitalStability(A, B, C, wn)
        isStable = true;
        Kn = ComputeKn(A, B, C, wn);
        if verbose
            fprintf('SUCCESS: b = %.2f mV (Stable, kn = %.2f, w = %.2f, f = %.2f Hz)\n', b*1000, Kn, wn, wn/(2*pi));
        end

    else
        isStable = false;
        Kn = 0;
        if verbose
            fprintf('FAIL: b = %.2f V (Unstable, w = %.2f, f = %.2f Hz)\n', b, wn, wn/(2*pi));
        end

    end

end

