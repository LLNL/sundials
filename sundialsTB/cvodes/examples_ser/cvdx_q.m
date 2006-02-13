function [qd, flag, new_data] = cvdx_q(t, y, data)
%CVDX_Q - quadrature function for the CVADX example problem.
%
%   See also: cvadx, CVQuadRhsFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/01/06 18:59:49 $


qd = y(3);

flag = 0;
new_data = [];
