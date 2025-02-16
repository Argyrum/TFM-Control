function J = lprsmatrdel(A, B, C, tau, w)
% Calculation of a point of the LPRS for matrix-vector system description having time delay "tau", dx/dt = Ax+Bu(t-tau); y = Cx
%   w - frequency

n = size(A, 1);
AINV = inv(A);
I = eye(n);
if w == 0
   J = (-0.5+1i*0.25*pi)*C*AINV*B;
else
   t = 2.*pi/w;
   AEXP = expm(0.5*A*t);
   AEXP2 = expm(A*t);
   AEXP3 = expm(A*(0.5*t-tau));
   re_lprs = -0.5*C*(AINV+t*inv(I-AEXP2)*AEXP3)*B;
   im_lprs = 0.25*pi*C*inv(I+AEXP)*(I+AEXP-2*AEXP3)*AINV*B;
   J = re_lprs+1i*im_lprs;
end