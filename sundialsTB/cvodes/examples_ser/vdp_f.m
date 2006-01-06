function [yd, new_data] = vdp_f(t, y, data)
%VDP_F - RHS function for the VDP example problem.
%
%   See also: vdp, CVRhsFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$


mu = data.mu;

yd = [            y(2)
        mu*(1-y(1)^2)*y(2)-y(1) ];

new_data = [];