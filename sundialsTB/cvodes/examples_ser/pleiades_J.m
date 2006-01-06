function [J, new_data] = pleiades_J(t, y, fy, data)
%PLEIADES_J - Jacobian function for the PLEIADES example problem.
%
%   see also: pleiades, CVDenseJacFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$

neq = 28;

J = zeros(neq,neq);
for i = 1:14
  J(i,14+i)=1.0;
end
for i = 2:7
  mi=i;
  for j = 1:i-1
    mj = j;
    rij = (y(i)-y(j))^2+(y(i+7)-y(j+7))^2;
    rij32 = rij^(3/2);
    rij52 = rij^(5/2);
    fjh = (1.0-3.0*(y(j)-y(i))^2/rij)/rij32;
    J(i+14,j)   = mj*fjh;
    J(j+14,i)   = mi*fjh;
    fjh = (1.0-3.0*(y(j+7)-y(i+7))^2/rij)/rij32;
    J(i+21,j+7) = mj*fjh;
    J(j+21,i+7) = mi*fjh;
    fjh = -3.0*(y(j)-y(i))*(y(j+7)-y(i+7))/rij52;
    J(i+14,j+7) = mj*fjh;
    J(j+14,i+7) = mi*fjh;
    J(i+21,j)   = mj*fjh;
    J(j+21,i)   = mi*fjh;
  end
end
for i = 1:7
  sumxx = 0.0;
  sumxy = 0.0;
  sumyy = 0.0;
  for j = 1:7
    if j ~= i
      sumxx = sumxx + J(i+14,j);
      sumxy = sumxy + J(i+14,j+7);
      sumyy = sumyy + J(i+21,j+7);
    end
  end
  J(i+14,i)   = -sumxx;
  J(i+14,i+7) = -sumxy;
  J(i+21,i)   = -sumxy;
  J(i+21,i+7) = -sumyy;
end

new_data = [];
