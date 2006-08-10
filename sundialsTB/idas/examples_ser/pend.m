function pend

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

IDAMalloc(@pend_f,t0,y0,yp0,options);

[status, y0_mod, yp0_mod] = IDACalcIC(tf, 'FindAlgebraic');
fprintf('Corrected IC\n');
disp([y0_mod yp0_mod])

it = 1;
time(it) = t0;
sol_y(it,:) = y0_mod';
sol_yp(it,:) = yp0_mod';
[pc(it) vc(it)] = pend_constr(t0,y0_mod);

t = t0;
t_start = clock;
while t < tf
  [status,t,y,yp] = IDASolve(tf,'OneStep');
  it = it+1;
  time(it) = t;
  sol_y(it,:) = y';
  sol_yp(it,:) = yp';
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



function [res, flag, new_data] = pend_f(t,y,yp)
%  Residual function F of The Simple Pendulum given as an index 1 DAE.
%  Automatically created by SimplePend_DAE1_IDA.m
%  Copyright (C) July 2006, Hannes Gruschinski
%  Rose-Hulman Institute of Technology, Terre Haute, IN

res = [  -yp(1)+y(3)  
	  -yp(2)+y(4)  
	  -yp(3)-2*y(1)*y(5)-3/10*y(3)  
	  -yp(4)+981/100-2*y(2)*y(5)-3/10*y(4)  
	  -2*y(5) + y(3)^2 - 3/10*y(1)*y(3) + y(4)^2 + y(2)*(981/100-3/10*y(4))
	];

flag = 0;
new_data = [];




function [J, flag, new_data] = pend_J(t,y,yp, rr, cj)
%  J the system iteration matrix of The Simple Pendulum given as an index 1 DAE.
%  Automatically created by SimplePend_DAE1_IDA.m
%  Copyright (C) July 2006, Hannes Gruschinski
%  Rose-Hulman Institute of Technology, Terre Haute, IN

J = [  -cj 0 1 0 0  
	  0 -cj 0 1 0  
	  -2*y(5) 0 -cj-3/10 0 -2*y(1)  
	  0 -2*y(5) 0 -cj-3/10 -2*y(2)  
	  -3/10*y(3) 981/100-3/10*y(4) 2*y(3)-3/10*y(1) 2*y(4)-3/10*y(2) -2
	];

flag = 0;
new_data = [];




function [pc, vc] = pend_constr(t,y)
%  Position and velocity constraints
%

pc = y(1)^2 + y(2)^2 - 1.0;
vc = y(1)*y(3) + y(2)*y(4);