function mcvsVanDPol_dns()
%mcvsVanDPol_dns - CVODES example problem (serial, dense)
%   van der Pol problem.

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


data.mu = 100.0;

t0 = 0.0;
tf = 300.0;
y0 = [2.0;0.0];

options = CVodeSetOptions('UserData',data,...
                          'RelTol',1.e-8,...
                          'AbsTol',1e-6,...
                          'JacobianFn',@djacfn);

mondata.mode = 'both';
mondata.skip = 20;
options = CVodeSetOptions(options,'MonitorFn',@CVodeMonitor,'MonitorData',mondata);

CVodeInit(@rhsfn, 'BDF', 'Newton', t0, y0, options);

ntout = 50;
dt = (tf-t0)/ntout;
tt = linspace(t0+dt,tf,ntout-1);

[status,t,y] = CVode(tt,'Normal');

CVodeFree;

figure;
plot(t,y(1,:),'r',t,y(1,:),'.');

% ===========================================================================

function [yd, flag, new_data] = rhsfn(t, y, data)
% Right-hand side function

mu = data.mu;

yd = [            y(2)
        mu*(1-y(1)^2)*y(2)-y(1) ];

flag = 0;
new_data = [];

return

% ===========================================================================

function [J, flag, new_data] = djacfn(t, y, fy, data)
% Dense Jacobian function (if using Newton)

mu = data.mu;

J = [         0                  1
      -2*mu*y(1)*y(2)-1    mu*(1-y(1)^2) ];

flag = 0;
new_data = [];