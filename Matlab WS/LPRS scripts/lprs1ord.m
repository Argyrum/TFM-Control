function J = lprs1ord(k, t, w)
% Calculation of a point of the LPRS for Transfer Function G(s) = k/(t*s+1),
%   w - frequency

if w == 0
   J = k*(0.5-1i*pi/4);
else
   al = pi/t/w;
   J = 0.5*k*(1-al*csch(al)-1i*0.5*pi*tanh(al/2));
end   