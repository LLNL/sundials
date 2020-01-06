function midasPendI1_dns
%midasPendI1_dns - Simple pendulum modeled as an index-1 DAE
%  The pendulum is modeled using the x and y positions with
%  the constraint x^2 + y^2 = L^2
%  The index-1 DAE formulation (in first-order form) includes
%  differential equations for the positions and velocities and
%  the acceleration-level constraint.

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


% x, y, vx, vy, lam
neq = 5;

t0 = 0.0;
tf = 10.0;

id = ones(neq,1);
id(5) = 0;

options = IDASetOptions('RelTol',1.e-6,...
                        'AbsTol',1.e-6,...
                        'VariableTypes',id,...
                        'MaxNumSteps', 1000,...
                        'LinearSolver','Dense',...
                        'JacobianFn',@pend_J);
%mondata.update = 100;
%options = IDASetOptions(options,'MonitorFn',@IDAMonitor,'MonitorData',mondata);

y0 = zeros(neq,1);
y0(1) = 1.0;
y0(5) = 0.1;
yp0 = zeros(neq,1);
fprintf('Estimated IC\n');
disp([y0 yp0])

IDAInit(@pend_f,t0,y0,yp0,options);

[status, y0_mod, yp0_mod] = IDACalcIC(tf, 'FindAlgebraic');
fprintf('Corrected IC\n');
disp([y0_mod yp0_mod])

it = 1;
time(it) = t0;
sol_y(it,:) = y0_mod';
[pc(it) vc(it)] = pend_constr(t0,y0_mod);

%t = t0;
%t_start = clock;
%while t < tf
%  [status,t,y] = IDASolve(tf,'OneStep');
%  it = it+1;
%  time(it) = t;
%  sol_y(it,:) = y';
%  % Compute position and velocity constraint violations
%  [pc(it) vc(it)] = pend_constr(t,y);
%end
%runtime = etime(clock,t_start);


dt = 0.1;
nt = ceil((tf-t0)/dt);

t_start = clock;
for it = 1:nt
  tout = t0 + it*dt;
  [status,t,y] = IDASolve(tout,'Normal');
  time(it) = t;
  sol_y(it,:) = y';
% Compute position and velocity constraint violations
  [pc(it) vc(it)] = pend_constr(t,y);
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
plot(time,sol_y(:,5));
box on
set(gca,'XLim',[t0 tf])
title('Lagrange multiplier');

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

% ================================================================================

function [res, flag, new_data] = pend_f(t,y,yp)
% Residual function for a simple pendulum
%   mass = 1.0
%   length = 1.0
%   damping coeff. = 0.3
%   g = 9.81

res = [
    -yp(1) + y(3)  
    -yp(2) + y(4)  
    -yp(3) - 2*y(1)*y(5) - 0.3*y(3)  
    -yp(4) + 9.81 - 2*y(2)*y(5) - 0.3*y(4)  
    -2*y(5) + y(3)^2 - 0.3*y(1)*y(3) + y(4)^2 + y(2)*(9.81-0.3*y(4))
	];

flag = 0;
new_data = [];

% ================================================================================

function [J, flag, new_data] = pend_J(t,y,yp, rr, cj)

J = [
    -cj 0 1 0 0  
    0 -cj 0 1 0  
    -2*y(5) 0 -cj-0.3 0 -2*y(1)  
    0 -2*y(5) 0 -cj-0.3 -2*y(2)  
    -0.3*y(3) 9.81-0.3*y(4) 2*y(3)-0.3*y(1) 2*y(4)-0.3*y(2) -2
    ];

flag = 0;
new_data = [];

% ================================================================================

function [pc, vc] = pend_constr(t,y)
%  Position and velocity constraints
%

pc = y(1)^2 + y(2)^2 - 1.0;
vc = y(1)*y(3) + y(2)*y(4);