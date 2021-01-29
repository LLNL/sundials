function mcvsRoberts_dns()
%mcvsRoberts_dns - CVODES example problem (serial, dense)
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
%   tolerance and a vector absolute tolerance and also illustrates
%   the rootfinding capability in CVODES.
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

data.p = [0.04; 1.0e4; 3.0e7];

t0 = 0.0;
y0 = [1.0;0.0;0.0];

options = CVodeSetOptions('UserData', data,...
                          'RelTol',1.e-8,...
                          'AbsTol',[1.e-8; 1.e-14; 1.e-6],...
                          'LinearSolver','Dense',...
                          'JacobianFn',@djacfn,...
                          'RootsFn',@rootfn, 'NumRoots',2);

mondata.sol = true;
mondata.mode = 'text';
mondata.skip = 9;
mondata.updt = 100;
options = CVodeSetOptions(options,'MonitorFn',@CVodeMonitor,'MonitorData',mondata);

CVodeInit(@rhsfn, 'BDF', 'Newton', t0, y0, options);

t1 = 0.4;
tmult = 10.0;
nout = 12;

iout = 0;
tout = t1;
while iout < nout

  [status,t,y] = CVode(tout,'Normal');
  
% Extract statistics
  si = CVodeGetStats;

% Print output
  fprintf('t = %0.2e   order = %1d  step = %0.2e',t, si.qlast, si.hlast);
  if(status == 2)
    fprintf(' ... Root found  %d   %d\n',si.RootInfo.roots(1), si.RootInfo.roots(2));
  else
    fprintf('\n');
  end
  fprintf('solution = [ %14.6e  %14.6e  %14.6e ]\n\n', y(1), y(2), y(3));

% Update output time
  if(status == 0)
    iout = iout+1;
    tout = tout*tmult;
  end
  
end

si = CVodeGetStats;

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

function [g, flag, new_data] = rootfn(t,y,data)
% Root finding function

g(1) = y(1) - 0.0001;
g(2) = y(3) - 0.01;

flag = 0;
new_data = [];
