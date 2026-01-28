function J = lprs2ord4(k, xi, w)
% Calculation of a point of the LPRS for Transfer Function G(s) = k*s/(s*s+2*xi*s+1),
%   w - frequency
%   xi > 1

if w == 0
   J = 0-1i*0;
else
   sq = sqrt(xi*xi-1);
   k1 = -0.5/sq;
   k2 = -k1;
   t1 = xi+sq;
   t2 = xi-sq;
   J = k*(lprs1ord(k1, t1, w)+lprs1ord(k2, t2, w));
end