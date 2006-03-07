function [yd, flag, new_data] = pvkx_fl(t, y, data)
%PVKX_FL - local RHS computation for the PVKX example problem.
%
%   See also: pvkx, CVGlocalFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.3 $Date: 2006/02/13 23:01:25 $

xmin  = data.xmin;
ml    = data.ml;
start = data.start;
dx    = data.dx;
Dc    = data.Dc;
yext  = data.yext;

for i = 2:ml(1)+1
  for j = 2:ml(2)+1
    for k = 2:ml(3)+1

      x = xmin + (start + [i-2 j-2 k-2] ) .* dx;
      v = velocity(x, data);
      s = source(x,data);

      [c, cl, cr] = stencil(yext,i,j,k);

      adv = v .* (cr-cl) ./ (2.0*dx);
      dif = Dc * (cr - 2.0*c + cl) / dx.^2;
    
      yd(i-1,j-1,k-1) = s + sum(dif-adv);
      
    end
  end
end

yd = reshape(yd,prod(ml),1);

flag = 0;
new_data = [];


function [c,cl,cr] = stencil(yext,i,j,k)

c = yext(i,j,k) * ones(1,3);
cl(1) = yext(i-1,j,  k  ); cr(1) = yext(i+1,j,  k  );
cl(2) = yext(i,  j-1,k  ); cr(2) = yext(i,  j+1,k  );
cl(3) = yext(i,  j,  k-1); cr(3) = yext(i,  j,  k+1); 


function v = velocity(x, data)

L = data.L;
Vc = data.Vc;
xmin = data.xmin;

y = x(2) - xmin(2) - L;

v(1) = Vc * (L+y) * (L-y);
v(2) = 0.0;
v(3) = 0.0;


function s = source(x, data)

A1 = data.A1;  A2 = data.A2;
S1 = data.S1;  S2 = data.S2;
X1 = data.X1;  X2 = data.X2;
Y1 = data.Y1;  Y2 = data.Y2;
Z1 = data.Z1;  Z2 = data.Z2;

s1 = A1 * exp(-(X1-x(1))^2/S1) * exp(-(Y1-x(2))^2/S1) * exp(-(Z1-x(3))^2/S1);
s2 = A2 * exp(-(X2-x(1))^2/S2) * exp(-(Y2-x(2))^2/S2) * exp(-(Z2-x(3))^2/S2);

s = s1 + s2;

if s < data.GMIN
  s = 0.0;
end
