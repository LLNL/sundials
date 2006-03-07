function [qd, flag, new_data] = cvbx_q(t, u, data)
%CVBX_Q - quadrature function for the CVBX example problem.
%
%   See also: cvbx, CVQuadRhsFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.3 $Date: 2006/02/13 23:01:27 $


mx = data.mx;
my = data.my;
dx = data.dx;
dy = data.dy;
xmax = data.xmax;
ymax = data.ymax;

qd1 = 0.0;
for j = 1:my
  for i = 1:mx
    uij = u(j+(i-1)*my);
    if j == 1 | j == mx
      del_y = dy/2;
    else
      del_y = dy;
    end
    if i == 1 | i == mx
      del_x = dx/2;
    else
      del_x = dx;
    end
    qd1 = qd1 + uij * del_x*del_y;
  end
end

qd1 = qd1 / (xmax*ymax);

qd(1) = qd1;

flag = 0;
new_data = [];
