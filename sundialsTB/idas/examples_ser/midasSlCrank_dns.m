function midasSlCrank_dns
%midasSlCrank_dns  Slider-crank example
% The multibody system consists of two bodies (crank and 
% connecting rod) with a translational-spring-damper (TSDA)
% and a constant force acting on the connecting rod. 
%
% The system has a single degree of freedom. It is modeled with 3
% generalized coordinates (crank angle, horizontal position of the
% translational joint, and angle of the connecting rod) and 
% therefore has 2 constraints.
%
% Example 6.1.8, pp. 271 in 
% Ed. Haug - Intermediate Dynamics, Prentiss Hall, 1992
%
% For its solution with IDAS, the resulting index-3 DAE is reformulated
% as a stabilized index-2 DAE (Gear-Gupta-Leimkhuler formulation) by 
% introducing 2 additional Lagrange multipliers and appending the
% velocity constraints.
%
%  |                                  |
%  |       /\                         |       /\
%  |      /  \                        |      /  \
%  | a/2 /    \ 1/2                   |     /    \
%  |    /      \                      |    /      \
%  |   /--TSDA--\                     |   /        \
%  |  /          \                    |  /          \
%  | / a/2        \ 1/2               | /            \---
%  |/              \                  |/  \ q1   q3 / \  \
%  +-----------------------------     +----------------------------
%                    \                |             \  |\
%                     \               |               -|-\
%                      \ 1            |                |  \
%                       \             |<----- q2 ----->|   \
%                        \                                  \
%                         \                                  \
%                          \ --> F                            \
%
% The local reference frame on the crank is positioned at the
% revolute joint on the ground. The crank has length a, mass m1, and
% intertia (with respect to the local frame) J1.
% The local reference frame on the conncting rod is positioned at the
% translational joint. The connecting rod has length 2, mass m2, and
% inertia J2.
% The TSDA has spring constant k, damping constant c, and free length l0.
% A constant horizontal force F acts on the connecting rod.

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
% $Revision$Date: 2007/10/26 16:30:48 $

% Problem data
data.a = 0.5;
data.J1 = 1.0;
data.m2 = 1.0;
data.J2 = 2.0;
data.k = 1.0;
data.c = 1.0;
data.l0 = 1.0;
data.F = 1.0;

% Integration limits
t0 = 0.0;
tf = 10.0;

% Specify algebraic variables
id = ones(10,1);
id(7:10) = 0.0;

% Integration options
options = IDASetOptions('UserData',data,...
                        'RelTol',1.e-6,...
                        'AbsTol',1.e-6,...
                        'VariableTypes',id,...
                        'suppressAlgVars','on',...
                        'MaxNumSteps', 1000,...
                        'LinearSolver','Dense');
% Set consistent IC
[yy0, yp0] = setIC(data);

% Initialize IDAS
IDAInit(@scRes,t0,yy0,yp0,options);

% Store inital time and IC
it = 1;
time(it) = t0;
pos(it,:) = yy0(1:3)';
vel(it,:) = yy0(4:6)';
acc(it,:) = yp0(4:6)';
lam(it,:) = yy0(7:8)';

% Compute constraint joint forces at initial time
Fc(it,:) = joint_forces(yy0, yp0, data);

% Call solver in ONE_STEP mode
t = t0;
while t < tf
  [status,t,yy] = IDASolve(tf,'OneStep');
  it = it+1;
  time(it) = t;
  pos(it,:) = yy(1:3)';
  vel(it,:) = yy(4:6)';
  lam(it,:) = yy(7:8)';
  yp = IDAGet('DerivSolution',t,1);
  acc(it,:) = yp(4:6)';
  Fc(it,:) = joint_forces(yy, yp, data);
end

fprintf('Solver stats:\n');
disp(IDAGetStats);

% Plot solution
figure;
set(gcf,'position',[475 250 1000 800]);

subplot(2,2,1)
hold on
plot(time,pos(:,1),'b');
plot(time,pos(:,2),'r');
plot(time,pos(:,3),'g');
box on
set(gca,'XLim',[t0 tf])
title('position');
legend('q','x','p');

subplot(2,2,2)
hold on
plot(time,vel(:,1),'b');
plot(time,vel(:,2),'r');
plot(time,vel(:,3),'g');
box on
set(gca,'XLim',[t0 tf])
title('velocity');
legend('q''', 'x''', 'p''');

subplot(2,2,3)
hold on
plot(time,acc(:,1),'b');
plot(time,acc(:,2),'r');
plot(time,acc(:,3),'g');
box on
set(gca,'XLim',[t0 tf])
title('acceleration');
legend('q''''', 'x''''', 'p''''');


subplot(2,2,4)
hold on
plot(time,lam(:,1),'b');
plot(time,lam(:,2),'r');
box on
set(gca,'XLim',[t0 tf])
title('Lagrange multipliers (cnstr. forces)');
legend('\lambda_1', '\lambda_2');


% Plot joint forces
figure;
set(gcf,'position',[275 150 800 800]);

plot(time,Fc);
set(gca,'XLim',[t0 tf])
title('joint forces');
legend('vrt. force in rev. crank-ground',...
       'hrz. force in rev. crank-ground',...
       'vrt. force in rev. crank-rod',...
       'hrz. force in rev. crank-rod',...
       'vrt. force in transl.',...
       'torque in transl.');


IDAFree;

% ====================================================================================
% Consistent IC
% ====================================================================================

function [yy, yp] = setIC(data)

a = data.a;
J1 = data.J1;
m2 = data.m2;
J2 = data.J2;

q = pi/2.0;    
p = -asin(a*sin(q));
x = cos(p) + a*cos(q);

yy = zeros(10,1);
yp = zeros(10,1);

yy(1) = q;  % crank angle 
yy(2) = x;  % slider position
yy(3) = p;  % conecting rod angle

Q = appl_forces(yy, data);
g = gamma(yy,data);
G = jac(yy,data);

% Extended mass matrix in index-1 formulation

MM = zeros(5,5);
MM(1,1) = J1;
MM(2,2) = m2;
MM(3,3) = J2;
MM(4:5,1:3) = G;
MM(1:3,4:5) = G';

% Right-hand side in index-1 formulation

b = [Q;g];

% Solution of MM*x = b

acc = MM^(-1)*b;

yp(4) = acc(1);
yp(5) = acc(2);
yp(6) = acc(3);
%
yy(7) = acc(4);
yy(8) = acc(5);

return

% ====================================================================================
% Constraint Jacobian
% ====================================================================================

function G = jac(yy,data)

a = data.a;

q = yy(1);  % crank angle 
x = yy(2);  % slider position
p = yy(3);  % conecting rod angle

qd = yy(4); % crank angular velocity
xd = yy(5); % slider horizontal velocity
pd = yy(6); % conecting rod angular velocity

s1 = sin(q);
c1 = cos(q);
s2 = sin(p);
c2 = cos(p);

G(1,1) = a*s1;
G(1,2) = 1.0;
G(1,3) = s2;

G(2,1) = -a*c1;
G(2,2) = 0.0;
G(2,3) = -c2;

return

% ====================================================================================
% Right-hand side of acceleration constraint
% ====================================================================================

function g = gamma(yy, data)

a = data.a;

q = yy(1);  % crank angle 
x = yy(2);  % slider position
p = yy(3);  % conecting rod angle

qd = yy(4); % crank angular velocity
xd = yy(5); % slider horizontal velocity
pd = yy(6); % conecting rod angular velocity

s1 = sin(q);
c1 = cos(q);
s2 = sin(p);
c2 = cos(p);

g(1,1) = - (a*qd^2*c1+pd^2*c2);
g(2,1) = - (a*qd^2*s1 + pd^2*s2);

return

% ====================================================================================
% Generalized applied forces calculation
% ====================================================================================

function Q = appl_forces(yy, data)

a = data.a;
k = data.k;
c = data.c;
l0 = data.l0;
F = data.F;

q = yy(1);  % crank angle 
x = yy(2);  % slider position
p = yy(3);  % conecting rod angle

qd = yy(4); % crank angular velocity
xd = yy(5); % slider horizontal velocity
pd = yy(6); % conecting rod angular velocity

s1 = sin(q);
c1 = cos(q);
s2 = sin(p);
c2 = cos(p);
s21 = s2*c1 - c2*s1;
c21 = c2*c1 + s2*s1;

l2 = x^2 - x*(c2+a*c1) + (1.0 + a^2)/4.0 + a*c21/2.0;
l = sqrt(l2);
ld = 2.0*x*xd - xd*(c2+a*c1) + x*(s2*pd+a*s1*qd) - a*s21*(pd-qd)/2.0;
ld = ld / (2.0*l);

f = k*(l-l0) + c*ld;
fl = f/l;

Q(1,1) = - fl * a * (s21/2.0 + x*s1) / 2.0;
Q(2,1) = fl * (c2/2.0 - x + a*c1/2.0) + F;
Q(3,1) = - fl * (x*s2 - a*s21/2.0) / 2.0 - F*s2;

return

% ====================================================================================
% Residual function
% ====================================================================================

function [res, flag, new_data] = scRes(t,yy,yp,data)

a = data.a;
J1 = data.J1;
m2 = data.m2;
J2 = data.J2;

q = yy(1);  % crank angle 
x = yy(2);  % slider position
p = yy(3);  % conecting rod angle

qd = yy(4); % crank angular velocity
xd = yy(5); % slider horizontal velocity
pd = yy(6); % conecting rod angular velocity

lam1 = yy(7); % Lagrange multiplier (cnstr)
lam2 = yy(8); % Lagrange multiplier (cnstr)

mu1 = yy(9);  % Lagrange multiplier (GGL)
mu2 = yy(10); % Lagrange multiplier (GGL)

s1 = sin(q);
c1 = cos(q);
s2 = sin(p);
c2 = cos(p);

% Generalized forces
Q = appl_forces(yy, data);

% Velocities (GGL modified)
res(1) = yp(1) - qd + a*s1*mu1 - a*c1*mu2;
res(2) = yp(2) - xd + mu1;
res(3) = yp(3) - pd + s2*mu1 - c2*mu2; 

% Dynamical equations
res(4) = J1*yp(4) - Q(1) + a*s1*lam1 - a*c1*lam2;
res(5) = m2*yp(5) - Q(2) + lam1;
res(6) = J2*yp(6) - Q(3) + s2*lam1 - c2*lam2; 

% Position constraints
res(7) = x - c2 - a*c1;
res(8) = -s2 - a*s1;

% Velocity constraints
res(9) = a*s1*qd + xd + s2*pd;
res(10) = -a*c1*qd - c2*pd;

flag = 0;
new_data = [];

return

% ====================================================================================
% Joint constraint forces
% ====================================================================================

function Fc = joint_forces(yy, yp, data)
% Compute joint reaction forces for given positins, velocities, and
% accelerations. This is done by including the reaction forces and torques,
% considering the free body diagrams for the crank and connecting rod, and
% writing the dynamical equilibrium equations.
%

a = data.a;
k = data.k;
c = data.c;
l0 = data.l0;
F = data.F;

J1 = data.J1;
m2 = data.m2;
J2 = data.J2;

q = yy(1);  % crank angle 
x = yy(2);  % slider position
p = yy(3);  % conecting rod angle

qd = yy(4); % crank angular velocity
xd = yy(5); % slider horizontal velocity
pd = yy(6); % conecting rod angular velocity

qdd = yp(4); % crank angular acc.
xdd = yp(5); % slider horizontal acc.
pdd = yp(6); % connecting rod angular acc.

s1 = sin(q);
c1 = cos(q);
s2 = sin(p);
c2 = cos(p);
s21 = s2*c1 - c2*s1;
c21 = c2*c1 + s2*s1;

l2 = x^2 - x*(c2+a*c1) + (1.0 + a^2)/4.0 + a*c21/2.0;
l = sqrt(l2);
ld = 2.0*x*xd - xd*(c2+a*c1) + x*(s2*pd+a*s1*qd) - a*s21*(pd-qd)/2.0;
ld = ld / (2.0*l);

f = k*(l-l0) + c*ld;
fl = f/l;

% TSDA forces acting on crank and connecting rod

Q1A(1) = x-0.5*c2-a/2*c1;
Q1A(2) = -0.5*(s2+a*s1);
Q1A(3) = -0.5*(a*x*s1+s21);

Q2A(1) = -Q1A(1);
Q2A(2) = -Q1A(2);
Q2A(3) = -0.5*(x*s2+a*s21);

Q1A = fl*Q1A';
Q2A = fl*Q2A';

QA = [-Q1A;-Q2A];

% Force F acting on connecting rod
FA = [0;0;0;F;0;-F*s2];

% Dynamic forces
MA = [0;0;J1*qdd;m2*xdd;0;J2*pdd];

% Dynamic equilibrium equations:
%  MA = QA+FA + FC
% where the joint constraint forces are:
%   Fc(1) = vertical force in revolute join crank-ground
%   Fc(2) = horizontal force in revolute join crank-ground
%   Fc(3) = vertical force in revolute joint crank-connecting rod
%   Fc(4) = vertical force in revolute joint crank-connecting rod
%   Fc(5) = vertical force in translational joint
%   Fc(6) = reaction torque acting in translational joint
% and therefore
%   FC = A * Fc
% where
%   A = [
%       0  1   0     1    0  0
%       1  0   1     0    0  0
%       0  0  a*c1 -1*s1  0  0
%       0  0   0    -1    0  0
%       0  0  -1     0    1  0
%       0  0  -c2   -s2   0  1
%       ];

b = MA-QA-FA;

Fc(4) = -b(4);
Fc(2) = b(1) - Fc(4);
if abs(c1) > 0.1
  Fc(3) = ( b(3) + a*s1*Fc(4) ) / (a*c1);
end
Fc(1) = b(2) - Fc(3); 
Fc(5) = b(5) + Fc(3);
Fc(6) = b(6) + c2*Fc(3) + s2*Fc(4);

return

