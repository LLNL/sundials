function [yd, flag, new_data] = pvnx_f(t, y, data)
%PVNX_F - RHS functin for the PVNX exampel problem.
%
%   see also: pvnx, CVRhsFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/01/06 18:59:46 $


alpha  = data.alpha;
nlocal = data.nlocal;
mype   = data.mype;

for i = 1:nlocal
  yd(i) = -alpha * (mype*nlocal + i) * y(i);
end

flag = 0;
new_data = [];
