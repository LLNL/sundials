function [fy, new_data] = kindiag_sys(y, data)
%KINDIAG_SYS - System function for the  KINDIAG example problem.
%
%   See also: kindiag, KINsysFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$

neq = length(y);
for i = 1:neq
  fy(i) = y(i)^2 - i^2;
end
new_data = [];