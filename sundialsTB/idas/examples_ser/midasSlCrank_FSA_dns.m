function midasSlCrank_FSA_dns
%midasSlCrank_FSA_dns - FSA for the slider-crank example
%
% Sensitivities w.r.t. k and c are computed
%
% See also: midasSlCrank_dns

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


% Problem data
data.a = 0.5;
data.J1 = 1.0;
data.m2 = 1.0;
data.J2 = 2.0;
data.l0 = 1.0;
data.F = 1.0;
data.params(1) = 1.0; % spring constant
data.params(2) = 1.0; % damper constant

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
IDAInit(@scRES,t0,yy0,yp0,options);

% FSA options
Ns = 2;
options = IDASensSetOptions('method','Simultaneous',...
                            'ErrControl',true,...
                            'ParamField','params',...
                            'ParamList',[1 2]);

% Sensitivity IC
yyS0 = zeros(10, Ns);
ypS0 = zeros(10, Ns);

% Initialize FSA
IDASensInit(Ns, [], yyS0, ypS0, options);

% Compute consistent IC 

% Store inital time and IC
it = 1;
time(it,1) = t0;
solution(it,:) = yy0';
sensitivity1(it,:) = yyS0(:,1)';
sensitivity2(it,:) = yyS0(:,2)';

% Call solver in ONE_STEP mode
t = t0;
while t < tf
  [status,t,y,yS] = IDASolve(tf,'OneStep');
  it = it+1;
  time(it,1) = t;
  solution(it,:) = y';
  sensitivity1(it,:) = yS(:,1)';
  sensitivity2(it,:) = yS(:,2)';
end

fprintf('Solver stats:\n');
disp(IDAGetStats);

IDAFree;

% Plot slider position and its sensitivities
figure;
set(gcf,'position',[475 250 1000 400]);

hold on
X = [time ; flipud(time)];
Y1 = [solution(:,2) ; flipud(solution(:,2)+sensitivity1(:,2))];
Y2 = [solution(:,2) ; flipud(solution(:,2)+sensitivity2(:,2))];


hp1 = patch(X,Y1,'r');
hp2 = patch(X,Y2,'b');

%set(hp1,'EdgeColor','none','FaceAlpha',0.5);
%set(hp2,'EdgeColor','none','FaceAlpha',0.5);

set(hp1,'EdgeColor','none');
set(hp2,'EdgeColor','none');

hp = plot(time,solution(:,2),'k');
set(hp,'LineWidth',2);

set(gca,'XLim',[t0 tf]);

box on
grid on


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

Q = force(yy, data);

yp(4) = Q(1)/J1;  % crank angular acceleration
yp(5) = Q(2)/m2;  % slider horizontal acceleration
yp(6) = Q(3)/J2;  % connecting rod angular acceleration

return

% ====================================================================================
% Generalized force calculation
% ====================================================================================

function Q = force(yy, data)

a = data.a;
k = data.params(1);
c = data.params(2);
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

Q(1) = - fl * a * (s21/2.0 + x*s1) / 2.0;
Q(2) = fl * (c2/2.0 - x + a*c1/2.0) + F;
Q(3) = - fl * (x*s2 - a*s21/2.0) / 2.0 - F*s2;

return

% ====================================================================================
% Residual function
% ====================================================================================

function [res, flag, new_data] = scRES(t,yy,yp,data)

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
Q = force(yy, data);

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

