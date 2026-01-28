function [wn, isStable, Kn] = FindWnPrecise(A, B, C, b, c, wini, verbose, precision)
% FINDWNPRECISE [wn, isStable, Kn] - Find wn such imag{J(wn)}=-pi*b/(4*c) within wini
%   Inverse of Compute_b

    arguments 
        A
        B
        C
        b
        c
        wini
        verbose = true
        precision = 100*eps
    end

    %options = optimset('Display', 'iter', 'TolX', precision);
    options = optimset('TolX', precision);
    obj_fun = @(w) abs(imag(lprsmatr(A, B, C, w)) + pi*b/(4*c));
    
    [wn, ~, exitflag] = fminsearch(obj_fun, wini, options);

    if (exitflag ~= 1)
        fprintf('ERROR - Could not find solution\n');

    elseif TestOrbitalStability(A, B, C, wn)
        isStable = true;
        Kn = ComputeKn(A, B, C, wn);
        if verbose
            fprintf('SUCCESS: wn = %.2f (Stable, kn = %.2f, f = %.2f Hz)\n', wn, Kn, wn/(2*pi));
        end

    else
        isStable = false;
        Kn = 0;
        if verbose
            fprintf('FAIL: wn = %.2f (Unstable, f = %.2f Hz)\n', wn, wn/(2*pi));
        end

    end
end

