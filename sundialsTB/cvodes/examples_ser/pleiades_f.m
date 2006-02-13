function [yd, flag, new_data] = pleiades_f(t, y, data)
%PLEIADES_F - RHS function for the PLEIADES example problems.
%
%   See also: pleiades, CVRhsFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/01/06 18:59:49 $


for i = 1:7
  sumx = 0.0;
  sumy = 0.0;
  for j = 1:7
    mj = j;
    rij = (y(i)-y(j))^2 + (y(i+7)-y(j+7))^2;
    rij32 = rij^(3/2);
    if j ~= i
      sumx = sumx + mj*(y(j)-y(i))/rij32;
      sumy = sumy + mj*(y(j+7)-y(i+7))/rij32;
    end
  end
  yd(i+14) = sumx;
  yd(i+21) = sumy;
end
for i = 1:14
  yd(i) = y(i+14);
end

flag = 0;
new_data = [];