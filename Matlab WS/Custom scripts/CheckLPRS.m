function [wS, iwS, KnS, bS, wE, KnE] = CheckLPRS(G, name, warr, c, b, verbose, plots)
%CHECKLPRS [wS, iwS, KnS, bS, wE, KnE] - Check stability equilibrium
%   Detailed explanation goes here

    arguments
        G
        name
        warr
        c = 0
        b = 0
        verbose = true
        plots = true
    end

    %% Compute LPRS for base plant G
    % Compute J(w)
    Gss = ss(G);
    J = arrayfun(@(w) lprsmatr(Gss.A, Gss.B, Gss.C, w), warr); % Complex - LPRS of G(w)
    imJ = imag(J);
    
    %icabove = find(1 == conv(imJ < pi*b/(4*c), [1 -1], 'same')); % Index - Crossings above -> below (index of above)
    %cabove = size(icabove, 2); % Nº of crossings above
    
    %icbelow = find(1 == conv(imJg > -pi*b/(4*c), [1 -1], 'same')); % Index - Crossings below -> above (index of below)
    %cbelow = size(icbelow, 2); % Nº of crossings below
    
    % Find the index(s) of the limit frequency(s)
    %JS = arrayfun(@(w) TestOrbitalStability(Gss.A, Gss.B, Gss.C, w), w); % Bool - Is J(w) stable
    %iwS = find(1 == conv(JS, [1 -1])) % Find transition 0 -> 1 (index of 0)
    
    % Find the index of the minimum limit stable frequency (first value from 0)
    [wS, iwS] = FindOrbitalStabilityLimit(Gss.A, Gss.B, Gss.C, warr); % First w to be stable

    if verbose, fprintf('\nLimit frequency for stability in %s:\n', name); end
    [bS, ~, KnS]= Compute_b(Gss.A, Gss.B, Gss.C, c, wS, verbose); % Find b and gain at stability limit

    if (c ~= 0) && (b ~= 0)
        if verbose, fprintf('\nEquilibrium frequency for b = %f, c = %i in %s:\n', b, c, name); end
        [wE, isStable, KnE] = FindWnPrecise(Gss.A, Gss.B, Gss.C, b, c, wS, verbose);
        if ~isStable
            wE = 0;
            KnE = 0;
        end
    else
        wE = 0;
        KnE = 0;
    end
    
    %% Plot LPRS information
    
    %disp('Rad/s of crossings at pi*b/(4*c):') <- Apparently, only crosses below
    %for wi = w(icabove)
        %fprintf('    w = %f (kn = %f, f = %.2f Hz)\n', wi, ComputeKn(Gss.A, Gss.B, Gss.C, wi), wi/(2*pi))
    %end
    %disp('Rad/s of crossings at -pi*b/(4*c):') <- Superseded by FindWnPrecise
    %for wi = w(icbelow)
        %fprintf('    w = %f (kn = %f, f = %.2f Hz)\n', wi, ComputeKn(Gss.A, Gss.B, Gss.C, wi), wi/(2*pi))
    %end

    if plots    
        figure;
        
        %Plot J
        subplot(2, 1, 1), plot(real(J), imJ);
        title('J(w) plot'), xlabel("Re J(w)"), ylabel("Im J(w)");
        if (c ~= 0) && (b ~= 0)
            yline(-pi*b/(4*c), 'r');
            legend(["J(w)", "-pi*b/(4*c)"]);
        end
        
        %Plot Im J vs w
        JS = 1:size(warr,2) >= iwS; % Stability mask for w (1 stable, 0 not)
        
        subplot(2, 1, 2), plot(warr(JS), imJ(JS), warr(~JS), imJ(~JS));
        title('Im J(w)'), xlabel("w"), ylabel("Im J(w)");
        if (c ~= 0) && (b ~= 0)
            yline(-pi*b/(4*c), 'r');
            legend(["Im J(w) Stable", "Im J(w) Unstable", "-pi*b/(4*c)"]);
        else
            legend(["Im J(w) Stable", "Im J(w) Unstable"]);
        end

        sgtitle(sprintf('LPRS characteristics for %s', name));
    end
end

