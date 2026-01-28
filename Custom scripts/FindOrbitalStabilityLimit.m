function [wS, iwS] = FindOrbitalStabilityLimit(A, B, C, wArr)
% FINDORBITALSTABILITYLIMIT [wS, iwS] - Search first stable w in wArr

    for i = 1:length(wArr)
        if TestOrbitalStability(A, B, C, wArr(i))
            iwS = i;
            wS = wArr(i);
            return
        end
    end
end

