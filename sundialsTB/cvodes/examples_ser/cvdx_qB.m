function [qBd, flag, new_data] = cvdx_qB(t, y, yB, data)
%CVDX_QB - adjoint quadrature function for the CVADX example problem.
%
%   See also: cavdx, CVQuadRhsFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/01/06 18:59:49 $

r1 = data.p(1);
r2 = data.p(2);
r3 = data.p(3);

y1 = y(1);
y2 = y(2);
y3 = y(3);

l1 = yB(1);
l2 = yB(2);
l3 = yB(3);

l21 = l2-l1;
l32 = l3-l2;
y23 = y2*y3;

qBd(1) = y1*l21;
qBd(2) = -y23*l21;
qBd(3) = l32*y2^2;

flag = 0;
new_data = [];