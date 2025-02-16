function J = lprs2ord3(k, w)
% Calculation of a point of the LPRS for Transfer Function G(s) = k*s/(s+1)^2,
%   w - frequency

if w == 0
   J = 0-1i*0;
else
   al = pi/w;
   chal = cosh(al);
   shal = sinh(al);
   J = k*(0.5*al*(-shal+al*chal)/shal/shal-1i*0.25*pi*al/(1+chal));
end