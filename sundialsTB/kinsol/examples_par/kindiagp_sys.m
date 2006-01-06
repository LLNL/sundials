function [fy, new_data] = kindiagp_sys(y, data)
%KINDIAGP_SYS - System function for the  KINDIAGP example problem.
%
%   See also: kindiagp, KINsysFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$

nlocal = data.nlocal;
mype = data.mype;
baseadd = mype * nlocal;

for i = 1:nlocal
  fy(i) = y(i)^2 - (i+baseadd)^2;
end
new_data = [];