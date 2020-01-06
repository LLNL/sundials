function midasPendI2_dns
%midasPendI1_dns - Simple pendulum modeled as an index-2 DAE
%  The pendulum is modeled using the x and y positions with
%  the constraint x^2 + y^2 = L^2
%  The stabilized index-2 (GGL formulation) DAE (in first-order form)
%  includes differential equations for the positions and velocities
%  with additional Lagrange multipliers included in the position
%  differential equations) and the position and velocity constraints.

% Radu Serban <radu@llnl.gov>
% SUNDIALS Copyright Start
% Copyright (c) 2002-2020, Lawrence Livermore National Security
% and Southern Methodist University.
% All rights reserved.
%
% See the top-level LICENSE and NOTICE files for details.
%
% SPDX-License-Identifier: BSD-3-Clause
% SUNDIALS Copyright End
% $Revision$Date: 2007/10/26 16:30:48 $

% x, y, vx, vy, lam, mu
neq = 6;

t0 = 0.0;
tf = 10.0;

id = ones(neq,1);
id(5) = 0;
id(6) = 0;

options = IDASetOptions('RelTol',1.e-6,...
                        'AbsTol',1.e-6,...
                        'VariableTypes',id,...
                        'suppressAlgVars','on',...
                        'MaxNumSteps', 1000,...
                        'LinearSolver','Dense');
y0 = zeros(neq,1);
yp0 = zeros(neq,1);
y0(1) = 1.0;
yp0(4) = 9.81;
fprintf('Consistent IC:\n');
disp([y0 yp0])

IDAInit(@pendGGL_f,t0,y0,yp0,options);

it = 1;
time(it) = t0;
sol_y(it,:) = y0';
[res, dummy1, status] = pendGGL_f(t0, y0, yp0);
pc(it) = res(5);
vc(it) = res(6);

t = t0;
t_start = clock;
while t < tf
  [status,t,y] = IDASolve(tf,'OneStep');
  it = it+1;
  time(it) = t;
  sol_y(it,:) = y';
  yp=yp0;
  % For verification purposes only, compute position and velocity constraint violations
  % (use dummy yp = yp0)
  [res, dummy1, status] = pendGGL_f(t, y, yp0);
  pc(it) = res(5);
  vc(it) = res(6);
  
end
runtime = etime(clock,t_start);

fprintf('Solver stats:\n');
disp(IDAGetStats);
fprintf('Run time: %f\n',runtime);


figure;

subplot(3,1,1)
hold on
plot(time,sol_y(:,1),'b');
plot(time,sol_y(:,2),'r');
box on
set(gca,'XLim',[t0 tf])
title('position');
legend('x','y');

subplot(3,1,2)
hold on
plot(time,sol_y(:,3),'b');
plot(time,sol_y(:,4),'r');
box on
set(gca,'XLim',[t0 tf])
title('velocity');
legend('v_x', 'v_y');

subplot(3,1,3)
hold on
plot(time,sol_y(:,5),'b');
plot(time,sol_y(:,6),'r');
box on
set(gca,'XLim',[t0 tf])
title('Lagrange multipliers');
legend('\lambda', '\mu');

figure

plotyy(time, pc, time, vc);
box on
title('position and velocity constraint violations');

figure

subplot(2,1,1)
plot(sol_y(:,1),sol_y(:,2));
axis equal
axis tight
box on
grid on
xlabel('x');
ylabel('y');
title('trajectory');

phi = atan2( sol_y(:,1) , sol_y(:,2) );
phi_d = ( sol_y(:,1).*sol_y(:,4) - sol_y(:,2).*sol_y(:,3) ) ./ ( sol_y(:,1).^2 + sol_y(:,2).^2 ) ;
subplot(2,1,2)
plot3(time,phi, phi_d);
xlabel('time');
ylabel('\phi');
zlabel('\phi^\prime');
view(-30,15);
set(gca,'XLim',[t0 tf])
grid on
box on
title('phase plot');

IDAFree;



function [res, flag, new_data] = pendGGL_f(t,yy,yp)

g = 9.81;
m = 1.0;
b = 0.3;
L = 1.0;

x = yy(1);    xd = yp(1);
y = yy(2);    yd = yp(2);
vx = yy(3);   vxd = yp(3);
vy = yy(4);   vyd = yp(4);

lam = yy(5);
mu = yy(6);

res(1) = -xd  + (vx+2*x*mu);
res(2) = -yd  + (vy+2*y*mu);
res(3) = -vxd + (-b*vx+2*x*lam)/m;
res(4) = -vyd + (m*g-b*vy+2*y*lam)/m;

res(5) = x^2 + y^2 - L^2;
res(6) = 2*x*vx + 2*y*vy;

flag = 0;
new_data = [];


