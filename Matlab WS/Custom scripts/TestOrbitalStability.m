function IsStable = TestOrbitalStability(A, B, C, w)
% Test whether system is locally orbitally asymptotically stable as described in LPRS, p26

    T = 2*pi/w; % Orbital period
    eath = exp(A*T/2); % Precompute e^(A*T/2) for reuse
    v = 2*(eye(size(eath))+eath)^(-1)*eath*B; % Velocity matrix for t = (T/2)-
    
    IsStable = false; 

    if C*v > 0 % Verify relay switch direction
        Phi_o = (eye(size(v*C))-(v*C)/(C*v))*eath; % Compute Phi_o matrix
        if ~sum( abs(eig(Phi_o)) >= 1 ) % Check all eigenvalues have mag < 1
            IsStable = true;
        end
    end
end

