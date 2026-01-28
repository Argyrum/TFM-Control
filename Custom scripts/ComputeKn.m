function Kn = ComputeKn(A, B, C, wn)
% COMPUTEKN Compute LPRS Relay Equivalent Gain (Kn) at wn

    Jw = lprsmatr(A, B, C, wn);
    Kn = -1/(2*real(Jw));

end

