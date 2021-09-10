function mcvsPleiades_non()
%mcvsPleiades_non - CVODES example problem (serial, nonstiff)

% Radu Serban <radu@llnl.gov>
% SUNDIALS Copyright Start
% Copyright (c) 2002-2021, Lawrence Livermore National Security
% and Southern Methodist University.
% All rights reserved.
%
% See the top-level LICENSE and NOTICE files for details.
%
% SPDX-License-Identifier: BSD-3-Clause
% SUNDIALS Copyright End
% $Revision$Date: 2007/08/21 23:09:18 $

neq = 28;

t0 = 0.0;
tout = 3.0;

y0 = zeros(neq,1);
y0(1)  =  3.0;
y0(2)  =  3.0;
y0(3)  = -1.0;
y0(4)  = -3.0;
y0(5)  =  2.0;
y0(6)  = -2.0;
y0(7)  =  2.0;
y0(8)  =  3.0;
y0(9)  = -3.0;
y0(10) =  2.0;
y0(13) = -4.0;
y0(14) =  4.0;
y0(20) =  1.75;
y0(21) = -1.5;
y0(25) = -1.25;
y0(26) =  1.0;

options = CVodeSetOptions('RelTol', 1.0e-7,...
                          'AbsTol', 1.0e-7,...
                          'StopTime',tout,...
                          'MaxNumSteps',2000);

CVodeInit(@rhsfn, 'Adams', 'Functional', t0, y0, options);

% Loop in one-step mode
t = t0;
i = 0;
while t < tout
  i = i+1;
  [status,t,y] = CVode(tout,'OneStep');
  time(i) = t;
  xx(:,i) = y(1:7);
  yy(:,i) = y(8:14);
end

% Display solver statistics
Stats = CVodeGetStats

% Free solver memory
CVodeFree;

% Plot body trajectories
colors = ['k','b','r','g','c','y','m'];
figure;
for i = 1:7
  plot(xx(i,:),yy(i,:),colors(i));
  hold on;
end
legend('Body 1','Body 2','Body 3','Body 4','Body 5','Body 6','Body 7');
title('Body Trajectories');
xlabel('x');
ylabel('y');
grid on;
axis square;

% ===========================================================================

function [yd, flag, new_data] = rhsfn(t, y, data)
% Right-hand side function

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

return

% ===========================================================================

function [J, flag, new_data] = djacfn(t, y, fy, data)
% Dense Jacobian function

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

flag = 0;
new_data = [];

return