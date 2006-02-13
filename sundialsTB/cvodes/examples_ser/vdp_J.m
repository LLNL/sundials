function [J, flag, new_data] = vdp_J(t, y, fy, data)
%VDP_J - Jacobian funciton for the VDP example problem
%
%   See also: vdp, CVDenseJacFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/01/06 18:59:49 $


mu = data.mu;

J = [         0                  1
      -2*mu*y(1)*y(2)-1    mu*(1-y(1)^2) ];

flag = 0;
new_data = [];