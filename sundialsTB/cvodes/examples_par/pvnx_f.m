function [yd, new_data] = pvnx_f(t, y, data)
%PVNX_F - RHS functin for the PVNX exampel problem.
%
%   see also: pvnx, CVRhsFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$


alpha  = data.alpha;
nlocal = data.nlocal;
mype   = data.mype;

for i = 1:nlocal
  yd(i) = -alpha * (mype*nlocal + i) * y(i);
end

new_data = [];
