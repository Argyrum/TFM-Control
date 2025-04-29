function CheckTF(TF, name, verbose, plots)
%CHECKTF Calculate plant information
    arguments
        TF
        name
        verbose = true
        plots = true
    end

    % Analyse plant
    if verbose
        fprintf('\nZeros in %s:\n', name);
        disp(zero(TF));
        fprintf('\nPoles in %s:\n', name);
        disp(pole(TF));
    end
    
    %% Plot plant
    if plots
        figure;
        subplot(2, 2, 1), step(TF);
        subplot(2, 2, 2), bode(TF);
        subplot(2, 2, 3), nyquist(TF);
        subplot(2, 2, 4), rlocusplot(TF);
        sgtitle(sprintf('%s response characteristics', name));
    end
end

