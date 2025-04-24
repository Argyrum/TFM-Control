function [wEP, isStable, Kn] = FindWnPrecise(A, B, C, b, c, wini, precision)
% Refines w at imag{J(w)}=-pi*b/(4*c) between wl(1) and wl(2) and computes
% kn after checking stability

    % Define precision as optional argument and check bounds were provided
    arguments A, B, C, b, c, wini,
        precision = 100*eps
    end

    %options = optimset('Display', 'iter', 'TolX', precision);
    options = optimset('TolX', precision);
    obj_fun = @(w) abs(imag(lprsmatr(A, B, C, w)) + pi*b/(4*c));
    
	%[wEP, error, exitflag] = fzero(@(w) imag(lprsmatr(A, B, C, w)) + pi*b/(4*c), wini, options);
    [wEP, ~, exitflag] = fminsearch(obj_fun, wini, options);

    if (exitflag ~= 1)
        fprintf('ERROR - Could not find solution\n');

    elseif TestOrbitalStability(A, B, C, wEP)
        isStable = true;
        Kn = ComputeKn(A, B, C, wEP);
        fprintf('SUCCESS: wn = %f (Stable, kn = %f, f = %.2f Hz)\n', wEP, Kn, wEP/(2*pi))

    else
        isStable = false;
        fprintf('FAIL: w = %f (Unstable, f = %.2f Hz)\n', wEP, wEP/(2*pi))

    end
end

