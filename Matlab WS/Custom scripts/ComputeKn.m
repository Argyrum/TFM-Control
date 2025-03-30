function [Kn] = ComputeKn(A, B, C, w)
% Compute LPRS Relay Equivalent Gain (Kn) for slow motion at w

Jw = lprsmatr(A, B, C, w);
Kn = -1/(2*real(Jw));

end

