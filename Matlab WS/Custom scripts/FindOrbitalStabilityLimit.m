function [wS, iwS] = FindOrbitalStabilityLimit(A, B, C, wArr)
% Searches and returns the first w in wArr that makes the SS orbitally stable.

    for i = 1:length(wArr)
        if TestOrbitalStability(A, B, C, wArr(i))
            iwS = i;
            wS = wArr(i);
            return
        end
    end
end

