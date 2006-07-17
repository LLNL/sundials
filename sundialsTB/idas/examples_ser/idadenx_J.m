function [J, flag, new_data] = idadenx_J(t, y, yp, rr, cj, data)
%IDADENX_J - Jacobian function for the IDADENX example problems.
%
%   see also: idadenx, IDADenseJacFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/02/13 23:01:27 $


r1 = data.p(1);
r2 = data.p(2);
r3 = data.p(3);

J(1,1) = -r1 - cj;
J(2,1) = r1;
J(3,1) = 1.0;

J(1,2) = r2*y(3);
J(2,2) = -r2*y(3) - 2*r3*y(2) - cj;
J(3,2) = 1.0;

J(1,3) = r2*y(2);
J(2,3) = -r2*y(2);
J(3,3) = 1.0;

flag = 0;
new_data = [];
