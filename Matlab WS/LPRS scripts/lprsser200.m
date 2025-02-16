function J = lprsser200(w, name, pr)
% Function calculating the LPRS at a given frequency based on the series formula (as a sum of 200 terms of the series)
%   'w' - current frequency,
%   'name' - name of m-file providing calculation of transfer function,
%   'pr' - parameters of transfer function

reloc = 0;
imloc = 0;
iodd = -1;

for k = 1:200
   iodd = -iodd;
   omk = k*w;
   reimloc = feval(name,omk,pr);
   reloc = reloc+iodd*real(reimloc);
   if iodd == 1
      imloc = imloc+imag(reimloc)/k;
   end
end
J = reloc+1i*imloc;
