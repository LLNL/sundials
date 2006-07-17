function [rr, flag, new_data] = idadenx_f(t, y, yp, data)
%IDADENX_F - Residual function for the IDADENX example problems.
%
%   See also: idadenx, IDAResFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/02/13 23:01:27 $


r1 = data.p(1);
r2 = data.p(2);
r3 = data.p(3);

rr(1) = -r1*y(1) + r2*y(2)*y(3) - yp(1);
rr(2) =  r1*y(1) - r2*y(2)*y(3) - r3*y(2)*y(2) - yp(2);
rr(3) = y(1) + y(2) + y(3) - 1.0;

flag = 0;
new_data = [];