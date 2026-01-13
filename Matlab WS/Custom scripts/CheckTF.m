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
        TFpoles = pole(TF);
        TFzeros = zero(TF);
        TFrd = length(TFpoles) - length(TFzeros);

        fprintf('\nZeros in %s:\n', name);
        disp(TFzeros);
        fprintf('\nPoles in %s:\n', name);
        disp(TFpoles);

        if max(real(TFpoles)) < 0
            fprintf("\nOK: %s is stable\n", name);
        else
            fprintf("\nERROR: %s is not stable\n", name);
        end

        if max(real(TFzeros)) < 0
            fprintf("\nOK: %s is minimum phase\n", name);
        else
            fprintf("\nWARN: %s is not minimum phase\n", name);
        end

        if TFrd > 1
            fprintf("\nWARN: %s is strict proper (rd > 1)\n", name);
        elseif TFrd == 1
            fprintf("\nOK: %s is strict proper (rd = 1)\n", name);
        elseif TFrd == 0
            fprintf("\nWARN: %s is proper (rd = 0)\n", name);
        else
            fprintf("\nERROR: %s is not causal (rd < 0)\n", name);
        end

        if TFrd == 1 && max(real(TFzeros)) < 0
            fprintf("\nOK: %s is ASPR\n", name);
        else
            fprintf("\nWARN: %s is not ASPR\n", name);
        end
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

