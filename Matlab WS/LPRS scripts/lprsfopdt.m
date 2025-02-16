function J = lprsfopdt(k, t, tau, w)
% Calculation of a point of the LPRS for transfer function G(s) = k*exp(-tau*s)/(t*s+1),
%   'w' - current frequency

if w == 0
   J = k*(0.5-1i*pi/4);
else
   al = pi/t/w;
   gm = tau/t;
   expal = exp(-al);
   expgm = exp(gm);
   J = 0.5*k*(1-al*expgm*csch(al)+1i*0.5*pi*(2*expal*expgm/(1+expal)-1));
end