function J = lprsmatrint(A, B, C, w)
% Calculation of a point of the LPRS for matrix-vector system description having an integrating plant, dx/dt = Ax+Bu; y = Cx
%   w - frequency

n = size(A, 1);
AINV = inv(A);
AINV2 = AINV*AINV;
I = eye(n);
if w == 0
   J = 0.5*C*AINV*B-1i*1000000.;
else
   t = 2.*pi/w;
   D = expm(0.5*A*t);
   re_lprs = 0.25*C*AINV2*(inv(I-D*D)*(D*D-(I+2.*t*A)*D+D*D*D-I)+D-I)*B;
   im_lprs = 0.0625*pi*C*AINV*B*t+0.125*pi*C*AINV*AINV*(inv(I-D*D)*(3*D*D-3*D-D*D*D+I)-D+I)*B;
   J = re_lprs+1i*im_lprs;
end