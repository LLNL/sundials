function [fy, flag, new_data] = kindiag_sys(y, data)
%KINDIAG_SYS - System function for the  KINDIAG example problem.
%
%   See also: kindiag, KINsysFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/01/06 19:00:06 $

neq = length(y);
for i = 1:neq
  fy(i) = y(i)^2 - i^2;
end
flag = 0;       % success
new_data = [];  % data was not modified