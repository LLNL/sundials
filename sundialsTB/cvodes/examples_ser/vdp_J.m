function [J, new_data] = vdp_J(t, y, fy, data)
%VDP_J - Jacobian funciton for the VDP example problem
%
%   See also: vdp, CVDenseJacFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$


mu = data.mu;

J = [         0                  1
      -2*mu*y(1)*y(2)-1    mu*(1-y(1)^2) ];

new_data = [];