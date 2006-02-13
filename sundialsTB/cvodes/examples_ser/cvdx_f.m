function [yd, flag, new_data] = cvdx_f(t, y, data)
%CVDX_F - RHS function for the CVDX, CVFDX, and CVADX example problems.
%
%   See also: cvdx, cvfdx, cvadx, CVRhsFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/01/06 18:59:49 $


r1 = data.p(1);
r2 = data.p(2);
r3 = data.p(3);

yd(1) = -r1*y(1) + r2*y(2)*y(3);
yd(3) = r3*y(2)*y(2);
yd(2) = -yd(1) - yd(3);

flag = 0;
new_data = [];