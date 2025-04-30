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

        if max(real(TFzeros)) < 0
            disp("OK: Is minimum phase");
        else
            disp("WARN: Is not minimum phase");
        end

        if max(real(TFpoles)) < 0
            disp("OK: Is stable");
        else
            disp("ERROR: Is not stable");
        end

        if TFrd > 1
            disp("WARN: Is strict proper (rd > 1)");
        elseif TFrd == 1
            disp("OK: Is strict proper (rd = 1)");
        elseif TFrd == 0
            disp("WARN: Is proper (rd = 0)");
        else
            disp("ERROR: Is not causal (rd < 0)");
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

