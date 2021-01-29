function mcvsHessian_FSA_ASA
%mcvsHessian_FSA_ASA - CVODES Hessian calculation example (FSA over ASA)
%  The following is a simple example illustrating the use
%  of simultaneous FSA and ASA computations in order to 
%  evaluate the Hessian of a functional depending on the
%  ODE solution. 
%
%  The forward problem consists of the following 3 equations
%
%    dy1/dt = - p1 * y1^2 - y3;
%    dy2/dt = - y2;
%    dy3/dt = - p2^2 * y2 * y3;
%
%  depending on the parameters p1 = 1.0 and p2 = 2.0.
%
%  The initial conditions are y1(0) = y2(0) = y390) = 1.0
%
%  The functional considered is
%
%           2
%          /
%   G(p) = |  0.5 * ( y1^2 + y2^2 + y3^2 ) dt
%          /
%          0
%  

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



Neq = 3;
Np  = 2;

t0 = 0.0;
tf = 2.0;

% ----------------------------------------
% User data structure
% ----------------------------------------

data.p = [1.0 2.0];

% ----------------------------------------
% Forward CVODES options
% ----------------------------------------

options = CVodeSetOptions('UserData',data,...
                          'RelTol',1.e-8,...
                          'AbsTol',1.0e-8,...
                          'LinearSolver','Dense');

optionsQ = CVodeQuadSetOptions('ErrControl',true,...
                               'RelTol',1.e-8,'AbsTol',1.e-8);

optionsS = CVodeSensSetOptions('method','Simultaneous',...
                               'ErrControl', true,...
                               'ParamScales', [1.0; 2.0]);

y0 = ones(Neq,1);
CVodeInit(@rhsfn, 'BDF', 'Newton', t0, y0, options);

q0 = 0.0;
CVodeQuadInit(@rhsQfn, q0, optionsQ);

s0 = zeros(Neq,Np);
CVodeSensInit(Np, @rhsSfn, s0, optionsS);

% ----------------------------------------
% Initialize ASA
% ----------------------------------------

CVodeAdjInit(100, 'Polynomial');

% ----------------------------------------
% Forward integration
% ----------------------------------------

fprintf('\nForward integration ');

[status, tret, y, q, yS] = CVode(tf,'Normal');
s = CVodeGetStats;

fprintf('(%d steps)\n', s.nst);
fprintf('G = %12.4e\n', q);

% ----------------------------------------
% Backward CVODES options
% ----------------------------------------

optionsB = CVodeSetOptions('UserData',data,...
                           'RelTol',1.e-8,...
                           'AbsTol',1.e-8,...
                           'LinearSolver','Dense',...
                           'SensDependent', true);

optionsQB = CVodeQuadSetOptions('ErrControl',true,...
                                'RelTol',1.e-8,'AbsTol',1.e-8,...
                                'SensDependent', true);


yB1 = zeros(2*Neq,1);
idxB1 = CVodeInitB(@rhsB1fn, 'BDF', 'Newton', tf, yB1, optionsB);

qB1 = zeros(2*Np,1);
CVodeQuadInitB(idxB1, @rhsQB1fn, qB1, optionsQB);


yB2 = zeros(2*Neq,1);
idxB2 = CVodeInitB(@rhsB2fn, 'BDF', 'Newton', tf, yB2, optionsB);

qB2 = zeros(2*Np,1);
CVodeQuadInitB(idxB2, @rhsQB2fn, qB2, optionsQB);

% ----------------------------------------
% Backward integration
% ----------------------------------------

fprintf('\nBackward integration ');

[status, tretB, yB, qB] = CVodeB(t0,'Normal');
sB1 = CVodeGetStatsB(idxB1);
sB2 = CVodeGetStatsB(idxB2);

fprintf('(%d steps pb1) (%d steps pb2)\n',sB1.nst, sB2.nst);

qB1 = -qB{idxB1};
qB2 = -qB{idxB2};

fprintf('Gradient\n');
fprintf('  %12.4e  %12.4e  (from backward pb. 1)\n',qB1(1:2));
fprintf('  %12.4e  %12.4e  (from backward pb. 2)\n',qB2(1:2));

fprintf('Hessian\n');
fprintf('  %12.4e  %12.4e\n',qB1(3), qB2(3));
fprintf('  %12.4e  %12.4e\n',qB1(4), qB2(4));

% ----------------------------------------
% Free memory
% ----------------------------------------

CVodeFree;

% ===========================================================================

function [yd, flag, new_data] = rhsfn(t, y, data)
% Right-hand side function

p1 = data.p(1);
p2 = data.p(2);

yd(1) = - p1 * y(1)^2 - y(3);
yd(2) = - y(2);
yd(3) = - p2^2 * y(2) * y(3);

flag = 0;
new_data = [];

return

% ===========================================================================

function [qd, flag, new_data] = rhsQfn(t, y, data)
% Forward quadrature integrand function

qd = 0.5 * ( y(1)^2 + y(2)^2 + y(3)^2 );

flag = 0;
new_data = [];

return

% ===========================================================================

function [ySd, flag, new_data] = rhsSfn(t,y,yd,yS,data)
% Sensitivity right-hand side function

p1 = data.p(1);
p2 = data.p(2);

s = yS(:,1);

fys1 = - 2.0*p1*y(1) * s(1) - s(3);
fys2 = - s(2);
fys3 = - p2*p2*y(3) * s(2) - p2*p2*y(2) * s(3);

ySd(1,1) = fys1 - y(1)*y(1);
ySd(2,1) = fys2;
ySd(3,1) = fys3;

s = yS(:,2);

fys1 = - 2.0*p1*y(1) * s(1) - s(3);
fys2 = - s(2);
fys3 = - p2*p2*y(3) * s(2) - p2*p2*y(2) * s(3);

ySd(1,2) = fys1;
ySd(2,2) = fys2;
ySd(3,2) = fys3 - 2.0*p2*y(2)*y(3);

flag = 0;
new_data = [];

return

% ===========================================================================

function [yBd, flag, new_data] = rhsB1fn(t, y, yS, yB, data)
% Backward problem right-hand side function for 1st Hessian column

p1 = data.p(1);
p2 = data.p(2);

y1 = y(1);
y2 = y(2);
y3 = y(3);

s1 = yS(1,1);
s2 = yS(2,1);
s3 = yS(3,1);

l1 = yB(1);
l2 = yB(2);
l3 = yB(3);

m1 = yB(4);
m2 = yB(5);
m3 = yB(6);

yBd(1) = 2.0*p1*y1 * l1     - y1;
yBd(2) = l2 + p2*p2*y3 * l3 - y2;
yBd(3) = l1 + p2*p2*y2 * l3 - y3;

yBd(4) = 2.0*p1*y1 * m1     + l1 * 2.0*(y1 + p1*s1) - s1;
yBd(5) = m2 + p2*p2*y3 * m3 + l3 * p2*p2*s3         - s2;
yBd(6) = m1 + p2*p2*y2 * m3 + l3 * p2*p2*s2         - s3;

flag = 0;
new_data = [];

return

% ===========================================================================

function [yBd, flag, new_data] = rhsB2fn(t, y, yS, yB, data)
% Backward problem right-hand side function 2nd Hessian column
 
p1 = data.p(1);
p2 = data.p(2);

y1 = y(1);
y2 = y(2);
y3 = y(3);

s1 = yS(1,2);
s2 = yS(2,2);
s3 = yS(3,2);

l1 = yB(1);
l2 = yB(2);
l3 = yB(3);

m1 = yB(4);
m2 = yB(5);
m3 = yB(6);

yBd(1) = 2.0*p1*y1 * l1     - y1;
yBd(2) = l2 + p2*p2*y3 * l3 - y2;
yBd(3) = l1 + p2*p2*y2 * l3 - y3;

yBd(4) = 2.0*p1*y1 * m1     + l1 * 2.0*p1*s1              - s1;
yBd(5) = m2 + p2*p2*y3 * m3 + l3 * (2.0*p2*y3 + p2*p2*s3) - s2;
yBd(6) = m1 + p2*p2*y2 * m3 + l3 * (2.0*p2*y3 + p2*p2*s2) - s3;

flag = 0;
new_data = [];

return

% ===========================================================================

function [qBd, flag, new_data] = rhsQB1fn(t, y, yS, yB, data)
% Backward problem quadrature integrand function for 1st Hessian column

p1 = data.p(1);
p2 = data.p(2);

y1 = y(1);
y2 = y(2);
y3 = y(3);

s1 = yS(1,1);
s2 = yS(2,1);
s3 = yS(3,1);

l1 = yB(1);
l2 = yB(2);
l3 = yB(3);

m1 = yB(4);
m2 = yB(5);
m3 = yB(6);

qBd(1) = -y1*y1 * l1;
qBd(2) = -2.0*p2*y2*y3 * l3;

qBd(3) = -y1*y1 * m1        - l1 * 2.0*y1*s1;
qBd(4) = -2.0*p2*y2*y3 * m3 - l3 * 2.0*(p2*y3*s2 + p2*y2*s3);

flag = 0;
new_data = [];

return

% ===========================================================================

function [qBd, flag, new_data] = rhsQB2fn(t, y, yS, yB, data)
% Backward problem quadrature integrand function for 2nd Hessian column

p1 = data.p(1);
p2 = data.p(2);

y1 = y(1);
y2 = y(2);
y3 = y(3);

s1 = yS(1,2);
s2 = yS(2,2);
s3 = yS(3,2);

l1 = yB(1);
l2 = yB(2);
l3 = yB(3);

m1 = yB(4);
m2 = yB(5);
m3 = yB(6);

qBd(1) = -y1*y1 * l1;
qBd(2) = -2.0*p2*y2*y3 * l3;

qBd(3) = -y1*y1 * m1        - l1 * 2.0*y1*s1;
qBd(4) = -2.0*p2*y2*y3 * m3 - l3 * 2.0*(p2*y3*s2 + p2*y2*s3 + y2*y3);

flag = 0;
new_data = [];

return