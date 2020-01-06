function midasRoberts_dns
%midasRoberts_dns - IDAS example problem (serial, dense)
%   The following is a simple example problem, with the coding
%   needed for its solution by IDAS. The problem is from
%   chemical kinetics, and consists of the following three rate
%   equations:         
%      dy1/dt = -.04*y1 + 1.e4*y2*y3
%      dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*(y2)^2
%           1 = y1 + y2 + y3
%   on the interval from t = 0.0 to t = 4.e10, with initial
%   conditions: y1 = 1.0, y2 = y3 = 0. The problem is stiff.
%   While integrating the system, we also use the rootfinding
%   feature to find the points at which y1 = 1e-4 or at which
%   y3 = 0.01.

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
% $Revision$Date: 2007/08/21 17:38:43 $

data.p = [0.04; 1.0e4; 3.0e7];

t0 = 0.0;
y0 = [1.0;0.0;0.0];
yp0 = [-0.04;0.04;0.0];

options = IDASetOptions('UserData', data,...
                        'RelTol',1.e-4,...
                        'AbsTol',[1.e-8; 1.e-14; 1.e-6],...
                        'LinearSolver','Dense',...
                        'JacobianFn',@djacfn);

options = IDASetOptions(options,'RootsFn',@rootfn, 'NumRoots',2);

%mondata.sol = true;
mondata.updt = 100;
options = IDASetOptions(options,'MonitorFn',@IDAMonitor,'MonitorData',mondata);

IDAInit(@resfn,t0,y0,yp0,options);

t1 = 0.4;
tmult = 10.0;
nout = 12;

fprintf('-----------------------------------------------------------------------\n');
fprintf('  t             y1           y2           y3');
fprintf('      | nst  k      h\n');
fprintf('-----------------------------------------------------------------------\n');

iout = 0;
tout = t1;
while iout < nout

  [status,t,y] = IDASolve(tout,'Normal');
  
% Extract statistics
  si = IDAGetStats;

% Print output
  if(status == 2)
    fprintf(' ... Root found  %d   %d\n',si.RootInfo.roots(1), si.RootInfo.roots(2));
  end
  fprintf('%10.4e %12.4e %12.4e %12.4e | %3d  %1d %12.4e\n',... 
          t, y(1), y(2), y(3), si.nst, si.qlast, si.hlast);

% Update output time
  if(status == 0)
    iout = iout+1;
    tout = tout*tmult;
  end
  
end

si = IDAGetStats;

IDAFree;

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

function [g, flag, new_data] = rootfn(t,y,yp,data)
% Root finding function

g(1) = y(1) - 0.0001;
g(2) = y(3) - 0.01;

flag = 0;
new_data = [];
