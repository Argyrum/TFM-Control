function J = lprsint(k, w)
% Calculation of a point of the LPRS for Transfer Function G(s) = k/s
%   w - frequency

if w == 0
   J = 0-1i*inf;
else   
   J = 0-1i*pi*pi*k/8/w;
end
