function [J, flag, new_data] = cvdx_J(t, y, fy, data)
%CVDX_J - Jacobian function for the CVDX, CVFDX, and CVADX example problems.
%
%   see also: cvdx, cvfdx, cvadx, CVDenseJacFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/01/06 18:59:48 $


r1 = data.p(1);
r2 = data.p(2);
r3 = data.p(3);

J(1,1) = -r1;
J(1,2) = r2*y(3);
J(1,3) = r2*y(2);

J(2,1) = r1;
J(2,2) = -r2*y(3) - 2*r3*y(2);
J(2,3) = -r2*y(2);

J(3,2) = 2*r3*y(2);

flag = 0;
new_data = [];
