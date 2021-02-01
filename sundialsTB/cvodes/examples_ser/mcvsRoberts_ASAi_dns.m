function mcvsRoberts_ASAi_dns()
%mcvsRoberts_ASAi_dns - CVODES adjoint sensitivity example problem (serial, dense)
%   The following is a simple example problem, with the coding
%   needed for its solution by CVODES. The problem is from
%   chemical kinetics, and consists of the following three rate
%   equations:         
%      dy1/dt = -p1*y1 + p2*y2*y3
%      dy2/dt =  p1*y1 - p2*y2*y3 - p3*(y2)^2
%      dy3/dt =  p3*(y2)^2
%   on the interval from t = 0.0 to t = 4.e10, with initial
%   conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff.
%
%   This program solves the problem with the BDF method,
%   Newton iteration with the CVDENSE dense linear solver, and a
%   user-supplied Jacobian routine. It uses a scalar relative 
%   tolerance and a vector absolute tolerance.
%
%   The gradient with respect to the problem parameters p1, p2,
%   and p3 of the following quantity:
%     G = int_t0^t1 y3(t) dt
%   is computed using ASA.
% 
%   Writing the original ODE as:
%     dy/dt = f(y,p)
%     y(t0) = y0(p)
%   then the gradient with respect to the parameters p of
%     G(p) = int_t0^t1 g(y,p) dt
%   is obtained as:
%     dG/dp = int_t0^t1 (g_p + lambda^T f_p ) dt + lambda^T(t0)*y0_p
%           = -xi^T(t0) + lambda^T(t0)*y0_p
%   where lambda and xi are solutions of:
%     d(lambda)/dt = - (f_y)^T * lambda - (g_y)^T
%     lambda(t1) = 0
%   and
%     d(xi)/dt = (g_p)^T + (f_p)^T * lambda
%     xi(t1) = 0
%   
%   During the forward integration, CVODES also evaluates G as
%     G = q(t1)
%   where
%     dq/dt = g(t,y,p)
%     q(t0) = 0

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


% ----------------------------------------
% User data structure
% ----------------------------------------

data.p = [0.04; 1.0e4; 3.0e7];

% ----------------------------------------
% Forward CVODES options
% ----------------------------------------


options = CVodeSetOptions('UserData',data,...
                          'RelTol',1.e-6,...
                          'AbsTol',[1.e-8; 1.e-14; 1.e-6],...
                          'LinearSolver','Dense',...
                          'JacobianFn',@djacfn);

mondata = struct;
mondata.sol = true;
mondata.mode = 'both';
options = CVodeSetOptions(options,...
                          'MonitorFn',@CVodeMonitor,...
                          'MonitorData',mondata);

t0 = 0.0;
y0 = [1.0;0.0;0.0];
CVodeInit(@rhsfn, 'BDF', 'Newton', t0, y0, options);


optionsQ = CVodeQuadSetOptions('ErrControl',true,...
                               'RelTol',1.e-6,'AbsTol',1.e-6);

q0 = 0.0;
CVodeQuadInit(@quadfn, q0, optionsQ);

% ----------------------------------------
% Initialize ASA
% ----------------------------------------

CVodeAdjInit(150, 'Hermite');

% ----------------------------------------
% Forward integration
% ----------------------------------------

fprintf('Forward integration ');

tout = 4.e7;
[status,t,y,q] = CVode(tout,'Normal');
s = CVodeGetStats;
fprintf('(%d steps)\n',s.nst);
fprintf('G = %12.4e\n',q);


fprintf('\nCheck point info\n');

ck = CVodeGet('CheckPointsInfo');
fprintf(['    t0         t1     nstep  order  step size\n']); 
for i = 1:length(ck)
  fprintf('%8.3e  %8.3e  %4d     %1d   %10.5e\n',...
          ck(i).t0, ck(i).t1, ck(i).nstep, ck(i).order, ck(i).step);
end
fprintf('\n');

% ----------------------------------------
% Backward CVODES options
% ----------------------------------------

optionsB = CVodeSetOptions('UserData',data,...
                           'RelTol',1.e-6,...
                           'AbsTol',1.e-8,...
                           'LinearSolver','Dense',...
                           'JacobianFn',@djacBfn);

mondataB = struct;
mondataB.mode = 'both';
optionsB = CVodeSetOptions(optionsB,...
                           'MonitorFn','CVodeMonitorB',...
                           'MonitorData', mondataB);

tB1 = 4.e7;
yB1 = [0.0;0.0;0.0];
idxB = CVodeInitB(@rhsBfn, 'BDF', 'Newton', tB1, yB1, optionsB);

optionsQB = CVodeQuadSetOptions('ErrControl',true,...
                                'RelTol',1.e-6,'AbsTol',1.e-3);

qB1 = [0.0;0.0;0.0];
CVodeQuadInitB(idxB, @quadBfn, qB1, optionsQB);

% ----------------------------------------
% Backward integration
% ----------------------------------------

fprintf('Backward integration ');

[status,t,yB,qB] = CVodeB(t0,'Normal');
sB=CVodeGetStatsB(idxB);
fprintf('(%d steps)\n',sB.nst);

fprintf('tB1:        %12.4e\n',tB1);
fprintf('dG/dp:      %12.4e  %12.4e  %12.4e\n',...
        -qB(1)+yB(1), -qB(2)+yB(2), -qB(3)+yB(3));
fprintf('lambda(t0): %12.4e  %12.4e  %12.4e\n',...
        yB(1),yB(2),yB(3));

% ----------------------------------------
% Free memory
% ----------------------------------------

CVodeFree;

% ===========================================================================

function [yd, flag, new_data] = rhsfn(t, y, data)
% Right-hand side function

r1 = data.p(1);
r2 = data.p(2);
r3 = data.p(3);

yd(1) = -r1*y(1) + r2*y(2)*y(3);
yd(3) = r3*y(2)*y(2);
yd(2) = -yd(1) - yd(3);

flag = 0;
new_data = [];

return

% ===========================================================================

function [qd, flag, new_data] = quadfn(t, y, data)
% Forward quadrature integrand function

qd = y(3);

flag = 0;
new_data = [];

return

% ===========================================================================

function [J, flag, new_data] = djacfn(t, y, fy, data)
% Dense Jacobian function

r1 = data.p(1);
r2 = data.p(2);
r3 = data.p(3);

J(1,1) = -r1;
J(1,2) = r2*y(3);
J(1,3) = r2*y(2);

J(2,1) = r1;
J(2,2) = -r2*y(3) - 2*r3*y(2);
J(2,3) = -r2*y(2);

J(3,2) = 2*r3*y(2);

flag = 0;
new_data = [];

return

% ===========================================================================

function [yBd, flag, new_data] = rhsBfn(t, y, yB, data)
% Backward problem right-hand side function

r1 = data.p(1);
r2 = data.p(2);
r3 = data.p(3);

y1 = y(1);
y2 = y(2);
y3 = y(3);

l1 = yB(1);
l2 = yB(2);
l3 = yB(3);

l21 = l2-l1;
l32 = l3-l2;
y23 = y2*y3;

yBd(1) = - r1*l21;
yBd(2) = r2*y3*l21 - 2.0*r3*y2*l32;
yBd(3) = r2*y2*l21 - 1.0;

flag = 0;
new_data = [];

return

% ===========================================================================

function [qBd, flag, new_data] = quadBfn(t, y, yB, data)
% Backward problem quadrature integrand function

r1 = data.p(1);
r2 = data.p(2);
r3 = data.p(3);

y1 = y(1);
y2 = y(2);
y3 = y(3);

l1 = yB(1);
l2 = yB(2);
l3 = yB(3);

l21 = l2-l1;
l32 = l3-l2;
y23 = y2*y3;

qBd(1) = y1*l21;
qBd(2) = -y23*l21;
qBd(3) = l32*y2^2;

flag = 0;
new_data = [];

return

% ===========================================================================

function [JB, flag, new_data] = djacBfn(t, y, yB, fyB, data)
% Backward problem Jacobian function

J = djacfn(t,y,[],data);
JB = -J';

flag = 0;
new_data = [];

return