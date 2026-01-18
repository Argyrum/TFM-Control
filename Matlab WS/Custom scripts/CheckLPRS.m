function [wS, iwS, KnS, bS, bE, KnE] = CheckLPRS(G, name, wArr, c, fe, verbose, plots)
%CHECKLPRS [wS, iwS, KnS, bS, bE, KnE] - Check stability and equilibrium
%   Detailed explanation goes here

    arguments
        G
        name
        wArr
        c = 0
        fe = 0
        verbose = true
        plots = true
    end

    %% Compute LPRS for base plant G
    % Compute J(w)
    Gss = ss(G);

    J = zeros(size(wArr));
    parfor idx = 1:size(wArr,2)
        J(idx) = lprsmatr(Gss.A, Gss.B, Gss.C, wArr(idx)); % Complex - LPRS of G(w)
    end

    imJ = imag(J);

    % Test stability and look for transitions
    JSMask = TestOrbitalStability(Gss.A, Gss.B, Gss.C, wArr);
    JSLimitsStartMask = find(conv(JSMask, [1 -1]) == 1);
    JSLimitsEndMask = find(conv(JSMask, [1 -1]) == -1);
    JSLimitsEndMask = JSLimitsEndMask(1:end-1) - 1; % Offset 1 forwards and discard last value to compensate for conv padding

    we = 2*pi*fe;

    if verbose
        fprintf('\nLimit frequencies for stability in %s:\n', name);
        for i=1:length(JSLimitsStartMask)
            fprintf('- Stability band %i -\n', i);
            Compute_b(Gss.A, Gss.B, Gss.C, c, wArr(JSLimitsStartMask(i)), true); % Find b and gain at stability limit
            if i <= length(JSLimitsEndMask)
                Compute_b(Gss.A, Gss.B, Gss.C, c, wArr(JSLimitsEndMask(i)), true); % Find b and gain at stability limit
            end
        end
    end

    % Ouput first stable frequency
    iwS = JSLimitsStartMask(1);
    wS = wArr(iwS);
    [bS, ~, KnS] = Compute_b(Gss.A, Gss.B, Gss.C, c, wS, false);

    if (c ~= 0) && (we ~= 0)
        if verbose, fprintf('\nHysteresis for equilibrium frequency at we = %f, c = %i in %s:\n', we, c, name); end
        [bE, isStable, KnE, JE] = Compute_b(Gss.A, Gss.B, Gss.C, c, we, verbose); % Find b and gain at objective we
        if ~isStable
            bE = 0;
            KnE = 0;
        end
    else
        bE = 0;
        KnE = 0;
    end
    
    %% Plot LPRS information

    if plots    
        figure;
        
        %Plot J
        subplot(2, 1, 1), hold on;
        plot(real(J), imJ, 'DisplayName', 'J(w)');
        title("J(w)"), xlabel("Re J(w)"), ylabel("Im J(w)");
        if (c ~= 0) && (bE ~= 0)
            yline(-pi*bE/(4*c), "--", "LineWidth", 2, 'DisplayName', '-pi*b/(4*c)');
            legend();
        end
        hold off;
        
        %Plot Im J vs w        
        subplot(2, 1, 2), hold on;
        plot(wArr(JSMask), imJ(JSMask), ".", 'DisplayName', 'Stable');
        plot(wArr(~JSMask), imJ(~JSMask), ".", 'DisplayName', 'Unstable');
        plot(wArr(JSLimitsStartMask), imJ(JSLimitsStartMask), "diamond", "LineWidth", 2, 'DisplayName', 'Stabilization');
        if ~isempty(JSLimitsEndMask)
            plot(wArr(JSLimitsEndMask), imJ(JSLimitsEndMask), "square", "LineWidth", 2, 'DisplayName', 'Destabilization');
        end
        if (c ~= 0) && (bE ~= 0)
            yline(-pi*bE/(4*c), ":", "LineWidth", 2, 'DisplayName', '-pi*bE/(4*c)');
            plot(we, imag(JE), "o", "LineWidth", 2, 'DisplayName', 'we')
        end

        title("Im J(w)"), xlabel("w"), ylabel("Im J(w)");
        legend();    
        
        xscale log;
        hold off;

        sgtitle(sprintf('LPRS characteristics for %s', name));
    end
end

