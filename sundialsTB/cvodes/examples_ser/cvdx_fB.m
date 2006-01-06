function [yBd, new_data] = cvdx_fB(t, y, yB, data)
%CVDX_FB - adjoint RHS functin for the CVADX example problem.
%
%   See alaos: cvadx, CVRhsFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$


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

yBd(1) = - r1*l21;
yBd(2) = r2*y3*l21 - 2.0*r3*y2*l32;
yBd(3) = r2*y2*l21 - 1.0;

new_data = [];