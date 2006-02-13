function [yd, flag, new_data] = vdp_f(t, y, data)
%VDP_F - RHS function for the VDP example problem.
%
%   See also: vdp, CVRhsFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/01/06 18:59:49 $


mu = data.mu;

yd = [            y(2)
        mu*(1-y(1)^2)*y(2)-y(1) ];

flag = 0;
new_data = [];