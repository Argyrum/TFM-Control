function J = lprs2ord2(k, xi, w)
% Calculation of a point of the LPRS for Transfer Function G(s) = k*s/(s*s+2*xi*s+1),
%   w - frequency
%   xi < 1

if w == 0
   J = 0-1i*0;
else
   al = pi*xi/w;
   sq = sqrt(1-xi*xi);
   bt = pi*sq/w;
   gm = al/bt;
   b = al*cos(bt)*sinh(al)+bt*sin(bt)*cosh(al);
   c = al*sin(bt)*cosh(al)-bt*cos(bt)*sinh(al);
   denom = sin(bt)^2+sinh(al)^2;
   J = 0.5*k*(-pi/w*sinh(al)*cos(bt)/denom+xi*(b+gm*c)/denom)-1i*0.25*k*pi/sq*sin(bt)/(cosh(al)+cos(bt));
end