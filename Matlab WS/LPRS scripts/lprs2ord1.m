function J = lprs2ord1(k, xi, w)
% Calculation of a point of the LPRS for Transfer Function G(s) = k/(s*s+2*xi*s+1),
%   w - frequency
%   xi < 1

if w == 0
   J = k*(0.5-1i*pi/4);
else
   al = pi*xi/w;
   sq = sqrt(1-xi*xi);
   bt = pi*sq/w;
   gm = al/bt;
   b = al*cos(bt)*sinh(al)+bt*sin(bt)*cosh(al);
   c = al*sin(bt)*cosh(al)-bt*cos(bt)*sinh(al);
   J = 0.5*k*(1-(b+gm*c)/(sin(bt)^2+sinh(al)^2))-1i*0.25*pi*k*(sinh(al)-gm*sin(bt))/(cosh(al)+cos(bt));
end