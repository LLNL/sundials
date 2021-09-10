function mcvsRoberts_FSA_dns()
%mcvsRoberts_FSA_dns - CVODES forward sensitivity example (serial, dense)
%   The following is a simple example problem, with the coding
%   needed for its solution by CVODES. The problem is from
%   chemical kinetics, and consists of the following three rate
%   equations:         
%      dy1/dt = -.04*y1 + 1.e4*y2*y3
%      dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*(y2)^2
%      dy3/dt = 3.e7*(y2)^2
%   on the interval from t = 0.0 to t = 4.e10, with initial
%   conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff.
%   While integrating the system, we also use the rootfinding
%   feature to find the points at which y1 = 1e-4 or at which
%   y3 = 0.01. This program solves the problem with the BDF method,
%   Newton iteration with the CVDENSE dense linear solver, and a
%   user-supplied Jacobian routine. It uses a scalar relative 
%   tolerance and a vector absolute tolerance.
%
%   Solution sensitivities with respect to the problem parameters
%   p1, p2, and p3 are computed using FSA. The sensitivity right-hand
%   side is given analytically through the user routine rhsSfn.
%   Tolerances for the sensitivity variables are estimated by
%   CVODES using the provided parameter scale information. The
%   sensitivity variables are included in the error test.
%
%   Output is printed in decades from t = .4 to t = 4.e10.
%   Run statistics (optional outputs) are printed at the end.

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

% -------------------
% User data structure
% -------------------

data.p = [0.04; 1.0e4; 3.0e7];

% ---------------------
% CVODES initialization
% ---------------------

options = CVodeSetOptions('UserData',data,...
                          'RelTol',1.e-4,...
                          'AbsTol',[1.e-8; 1.e-14; 1.e-6],...
                          'JacobianFn',@djacfn);

mondata = struct;
mondata.mode = 'both';
mondata.sol = true;
mondata.sensi = true;
options = CVodeSetOptions(options,'MonitorFn',@CVodeMonitor,'MonitorData',mondata);

t0 = 0.0;
y0 = [1.0;0.0;0.0];

CVodeInit(@rhsfn, 'BDF', 'Newton', t0, y0, options);


% ------------------
% FSA initialization
% ------------------

Ns = 2;
yS0 = zeros(3,Ns);

% Case 1: user-provided sensitivity RHS

FSAoptions = CVodeSensSetOptions('method','Simultaneous',...
                                 'ErrControl', true,...
                                 'ParamScales', [0.04; 1.0e4]);
CVodeSensInit(Ns, @rhsSfn, yS0, FSAoptions);

% Case 2: internal DQ approximation

%FSAoptions = CVodeSensSetOptions('method','Simultaneous',...
%                                 'ErrControl', true,...
%                                 'ParamField', 'p',...
%                                 'ParamList', [1 2],...
%                                 'ParamScales', [0.04 1.0e4]);
%CVodeSensInit(Ns, [], yS0, FSAoptions);

% ----------------
% Problem solution
% ----------------

t1 = 0.4;
tmult = 10.0;
nout = 12;

iout = 0;
tout = t1;
while 1,
  [status, t, y, yS] = CVode(tout,'Normal');
  fprintf('t = %0.2e\n',t);
  fprintf('solution      = [ %14.6e  %14.6e  %14.6e ]\n', y(1), y(2), y(3));
  fprintf('sensitivity 1 = [ %14.6e  %14.6e  %14.6e ]\n', yS(1,1), yS(2,1), yS(3,1));
  fprintf('sensitivity 2 = [ %14.6e  %14.6e  %14.6e ]\n\n', yS(1,2), yS(2,2), yS(3,2));
  if(status ==0)
    iout = iout+1;
    tout = tout*tmult;
  end
  if iout==nout
    break;
  end
end

si = CVodeGetStats

% -----------
% Free memory
% -----------

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

function [ySd, flag, new_data] = rhsSfn(t,y,yd,yS,data)
% Sensitivity right-hand side function

r1 = data.p(1);
r2 = data.p(2);
r3 = data.p(3);

% r1

yS1 = yS(:,1);
yS1d = zeros(3,1);

yS1d(1) = -r1*yS1(1) + r2*y(3)*yS1(2) + r2*y(2)*yS1(3);
yS1d(3) = 2*r3*y(2)*yS1(2);
yS1d(2) = -yS1d(1)-yS1d(3);

yS1d(1) = yS1d(1) - y(1);
yS1d(2) = yS1d(2) + y(1);

% r2

yS2 = yS(:,2);
yS2d = zeros(3,1);

yS2d(1) = -r1*yS2(1) + r2*y(3)*yS2(2) + r2*y(2)*yS2(3);
yS2d(3) = 2*r3*y(2)*yS2(2);
yS2d(2) = -yS2d(1)-yS2d(3);

yS2d(1) = yS2d(1) + y(2)*y(3);
yS2d(2) = yS2d(2) - y(2)*y(3);

% Return values

ySd = [yS1d yS2d];
flag = 0;
new_data = [];

return
