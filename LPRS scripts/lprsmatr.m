function J = lprsmatr(A, B, C, w)
% Calculation of a point of the LPRS for matrix-vector system description, dx/dt = Ax+Bu; y = Cx
%   w - frequency

n = size(A,1);
AINV = inv(A);
I = eye(n);
if w == 0
  J = (-0.5+1i*0.25*pi)*C*AINV*B;
else
  t = 2.*pi/w;
  AEXP = expm(0.5*A*t);
  AEXP2 = expm(A*t);
  re_lprs = -0.5*C*(AINV+t*inv(I-AEXP2)*AEXP)*B;
  im_lprs = 0.25*pi*C*inv(I+AEXP)*(I-AEXP)*AINV*B;
  J = re_lprs+1i*im_lprs;
end

% n = size(A,1);
% I = eye(n);
% if w == 0
%    J = (-0.5+1i*0.25*pi)*C/(A)*B;
% else
%    t = 2.*pi/w;
%    AEXP = expm(0.5*A*t);
%    AEXP2 = expm(A*t);
%    re_lprs = -0.5*C*( inv(A) + t*(I-AEXP2)\AEXP )*B;
%    im_lprs = 0.25*pi*C/(I+AEXP)*(I-AEXP)/(A)*B;
%    J = re_lprs + 1i*im_lprs;
% end