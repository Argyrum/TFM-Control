function [isStable, isMP, TFrd, isASPR] = CheckTF(TF, name, verbose, plots)
%CHECKTF Calculate plant information
    arguments
        TF
        name
        verbose = true
        plots = true
    end

    % Analyse plant
        TFpoles = pole(TF);
        isStable = max(real(TFpoles)) < 0;
        TFzeros = zero(TF);
        isMP = max(real(TFzeros)) < 0;
        TFrd = length(TFpoles) - length(TFzeros);
        isASPR = TFrd == 1 && isMP;

    if verbose
        fprintf('\nZeros in %s:\n', name);
        disp(TFzeros);
        fprintf('\nPoles in %s:\n', name);
        disp(TFpoles);

        if isStable
            fprintf("\nOK: %s is stable\n", name);
        else
            fprintf("\nERROR: %s is not stable\n", name);
        end

        if isMP
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

        if isASPR
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

