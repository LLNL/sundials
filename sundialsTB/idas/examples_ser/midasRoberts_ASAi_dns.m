function midasRoberts_ASAi_dns()
%midasRoberts_ASAi_dns - IDAS ASA example problem (serial, dense)
%   The following is a simple example problem, with the coding
%   needed for its solution by IDAS. The problem is from
%   chemical kinetics, and consists of the following three rate
%   equations:         
%      dy1/dt = -p1*y1 + p2*y2*y3
%      dy2/dt =  p1*y1 - p2*y2*y3 - p3*(y2)^2
%           0 = y1 + y2 + y3 - 1
%   on the interval from t = 0.0 to t = 4.e7, with initial
%   conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff.
%   While integrating the system, we also use the rootfinding
%   feature to find the points at which y1 = 1e-4 or at which
%   y3 = 0.01.
%
%   The gradient with respect to the problem parameters p1, p2,
%   and p3 of the following quantity:
%     G = int_t0^t1 y3(t) dt
%   is computed using ASA.
%          
%   The gradient dG/dp is obtained as:
%     dG/dp = [ int_t0^tf y1*(l1-l2) dt ,
%               int_t0^tf -y2*y3*(l1-l2) dt , 
%               int_t0^tf y2^2*l2 dt          ]
%              
%   where l = [l1, l2, l3] is solutions of:
%     dl1/dt = p1*l1 - p1*l2 + l3
%     dl2/dt = -p2*y3*l1 + (p2*y3+2*p3*y2)*l2 + l3
%          0 = -p2*y2*l1 + p2*y2*l2 + l3 + 1
%   with final conditions
%     l1(tf) = l2(tf) = 0.0  and l3(tf) = -1.0
%
%   All integrals (appearing in G and dG/dp) are computed using
%   the quadrature integration features in IDAS.

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
% $Revision$Date: 2007/08/21 17:38:43 $

% Problem parameters
% ------------------

data.p = [0.04; 1.0e4; 3.0e7];

% Initialize forward problem
% --------------------------

options = IDASetOptions('UserData', data,...
                        'RelTol',1.e-4,...
                        'AbsTol',[1.e-8; 1.e-14; 1.e-6],...
                        'LinearSolver','Dense',...
                        'JacobianFn',@djacfn);

%mondata.sol = true;
%mondata.updt = 100;
%options = IDASetOptions(options,'MonitorFn',@IDAMonitor,'MonitorData',mondata);

t0 = 0.0;
y = [1.0;0.0;0.0];
yp = [-0.04;0.04;0.0];
IDAInit(@resfn,t0,y,yp,options);

% Initialize forward quadrature (G)
% ---------------------------------

optionsQ = IDAQuadSetOptions('ErrControl',true,...
                             'RelTol',1.e-4,'AbsTol',1.e-6);
q = 0.0;
IDAQuadInit(@quadfn, q, optionsQ);

% Activate ASA
% ------------

IDAAdjInit(150, 'Hermite');

% Forward integration
% -------------------

fprintf('Forward integration ');
tf =  4.e7;
[status, t, y, q] = IDASolve(tf,'Normal');
si = IDAGetStats;
fprintf('(%d steps)\n',si.nst);

fprintf('G = %12.4e\n',q(1));

% Initialize backward problem
% ---------------------------

optionsB = IDASetOptions('UserData',data,...
                         'RelTol',1.e-6,...
                         'AbsTol',1.e-3,...
                         'LinearSolver','Dense');
%mondataB = struct;
%optionsB = IDASetOptions(optionsB,'MonitorFn',@IDAMonitorB,'MonitorData',mondataB);

yB = [0.0 ; 0.0 ; -1.0];
yBp = [ -1.0 ; -1.0 ; 0.0 ];
idxB = IDAInitB(@resfnB,tf,yB,yBp,optionsB);

% Initialize backward quadratures (dG/dp)
% ---------------------------------------

optionsQB = IDAQuadSetOptions('ErrControl',true,...
                              'RelTol',1.e-6,'AbsTol',1.e-3);

qB = [0.0;0.0;0.0];
IDAQuadInitB(idxB, @quadfnB, qB, optionsQB);

% Backward integration
% --------------------

fprintf('Backward integration ');
[status, t, yB, qB] = IDASolveB(t0,'Normal');
siB = IDAGetStatsB(idxB);
fprintf('(%d steps)\n',siB.nst);

fprintf('dG/dp:      %12.4e  %12.4e  %12.4e\n',...
        -qB(1),-qB(2),-qB(3));
fprintf('lambda(t0): %12.4e  %12.4e  %12.4e\n',...
        yB(1),yB(2),yB(3));

% Free IDAS memory
% ----------------

IDAFree;

return

% ===========================================================================
function [rr, flag, new_data] = resfn(t, y, yp, data)
% DAE residual function

r1 = data.p(1);
r2 = data.p(2);
r3 = data.p(3);

rr(1) = -r1*y(1) + r2*y(2)*y(3) - yp(1);
rr(2) =  r1*y(1) - r2*y(2)*y(3) - r3*y(2)*y(2) - yp(2);
rr(3) = y(1) + y(2) + y(3) - 1.0;

flag = 0;
new_data = [];

% ===========================================================================
function [J, flag, new_data] = djacfn(t, y, yp, rr, cj, data)
% Dense Jacobian function

r1 = data.p(1);
r2 = data.p(2);
r3 = data.p(3);

J(1,1) = -r1 - cj;
J(2,1) = r1;
J(3,1) = 1.0;

J(1,2) = r2*y(3);
J(2,2) = -r2*y(3) - 2*r3*y(2) - cj;
J(3,2) = 1.0;

J(1,3) = r2*y(2);
J(2,3) = -r2*y(2);
J(3,3) = 1.0;

flag = 0;
new_data = [];

% ===========================================================================
function [qd, flag, new_data] = quadfn(t, y, yp, data)
% Forward quadrature integrand function

qd = y(3);

flag = 0;
new_data = [];

return

% ===========================================================================
function [rrB, flag, new_data] = resfnB(t, y, yp, yB, yBp, data)
% Adjoint residual function

r1 = data.p(1);
r2 = data.p(2);
r3 = data.p(3);


rrB(1) = yBp(1) - r1*(yB(1)-yB(2)) - yB(3);
rrB(2) = yBp(2) + r2*y(3)*(yB(1)-yB(2)) - 2.0*r3*y(2)*yB(2) - yB(3);
rrB(3) = -r2*y(2)*(yB(1)-yB(2)) + yB(3) + 1.0;

flag = 0;
new_data = [];

return

% ===========================================================================
function [qBd, flag, new_data] = quadfnB(t, y, yp, yB, ypB, data)
% Backward problem quadrature integrand function

qBd(1) = y(1)*(yB(1)-yB(2));
qBd(2) = -y(2)*y(3)*(yB(1)-yB(2));
qBd(3) = y(2)^2*yB(2);

flag = 0;
new_data = [];

return

