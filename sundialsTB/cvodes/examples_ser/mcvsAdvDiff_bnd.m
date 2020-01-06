function mcvsAdvDiff_bnd()
%mcvsAdvDiff_bnd - CVODES example problem (serial, band)
%   The following is a simple example problem with a banded Jacobian,
%   with the program for its solution by CVODES.
%   The problem is the semi-discrete form of the advection-diffusion
%   equation in 2-D:
%      du/dt = d^2 u / dx^2 + .5 du/dx + d^2 u / dy^2
%   on the rectangle 0 <= x <= 2, 0 <= y <= 1, and the time
%   interval 0 <= t <= 1. Homogeneous Dirichlet boundary conditions
%   are posed, and the initial condition is
%      u(x,y,t=0) = x(2-x)y(1-y)exp(5xy).
%   The PDE is discretized on a uniform MX+2 by MY+2 grid with
%   central differencing, and with boundary values eliminated,
%   leaving an ODE system of size NEQ = MX*MY.
%   This program solves the problem with the BDF method, Newton
%   iteration with the CVBAND band linear solver, and a user-supplied
%   Jacobian routine.
%   It uses scalar relative and absolute tolerances.
%   Output is printed at t = .1, .2, ..., 1.
%   Run statistics (optional outputs) are printed at the end.

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
% $Revision$Date: 2007/08/21 23:09:18 $

xmax = 2.0;
ymax = 1.0;
mx = 10;
my = 5;

rtol = 0.0;
atol = 1.0e-5;

t0 = 0.0;
dtout = 0.1;
nout = 6;

dx = xmax/(mx+1);
dy =  ymax/(my+1);

% Problem data structure
data.xmax = xmax;
data.ymax = ymax;
data.mx = mx;
data.my = my;
data.dx = dx;
data.dy = dy;
data.hdcoef = 1.0/dx^2;
data.hacoef = 0.5/(2.0*dx);
data.vdcoef = 1.0/dy^2;

% Options for integration
options = CVodeSetOptions('UserData',data,...
                          'RelTol',rtol,...
                          'AbsTol',atol,...
                          'LinearSolver','Band',...
                          'JacobianFn',@bjacfn,...
                          'UpperBwidth',my,...
                          'LowerBwidth',my);

mondata.grph = false;
options = CVodeSetOptions(options,...
                          'MonitorFn',@CVodeMonitor,...
                          'MonitorData',mondata);

% Initial conditions for states
t = t0;
u = zeros(mx*my,1);
for j = 1:my
  y = j * data.dy;
  for i = 1:mx
    x = i * data.dx;
    u(j + (i-1)*my) = x*(xmax-x)*y*(ymax-y)*exp(5.0*x*y);
  end
end

% Initialize integrator
CVodeInit(@rhsfn, 'BDF', 'Newton', t, u, options);

% Initialize quadratures (with default optional inputs)
q = 0.0;
CVodeQuadInit(@quadfn, q);


ff=figure;
hold on;
box

umax = norm(u,'inf');
uavg = quadfn(t,u,data);
fprintf('At t = %f   max.norm(u) = %e\n',t, umax);

for i = 1:nout

  t_old = t;
  uavg_old = uavg;

  tout = t + dtout;
  [status,t,u, q] = CVode(tout, 'Normal');

  if status ~= 0
    return
  end
  
  uavg = quadfn(t,u,data);
  umax = norm(u,'inf');
  fprintf('At t = %f   max.norm(u) = %e\n',tout, umax);

  figure(ff);
  plot([t_old t],[uavg_old uavg]);
  plot([t0 tout], [q q]/(tout-t0), 'r-');
  plot([tout tout], [0 q]/(tout-t0), 'r-');
  
end
  
si= CVodeGetStats

CVodeFree;

return

% ===========================================================================

function [ud, flag, new_data] = rhsfn(t, u, data)
% Right-hand side function

mx = data.mx;
my = data.my;
hordc = data.hdcoef;
horac = data.hacoef;
verdc = data.vdcoef;

for j = 1:my
  for i = 1:mx
    uij = u(j+(i-1)*my);
    if j == 1
      udn = 0.0;
    else
      udn = u(j-1+(i-1)*my);
    end
    if j == my
      uup = 0.0;
    else
      uup = u(j+1+(i-1)*my);
    end
    if i == 1
      ult = 0.0;
    else
      ult = u(j+(i-2)*my);
    end
    if i == mx
      urt = 0.0;
    else
      urt = u(j+i*my);
    end
    
    hdiff = hordc * (ult - 2*uij + urt);
    hadv = horac * (urt - ult);
    vdiff = verdc * (uup - 2*uij + udn);
    ud(j+(i-1)*my) = hdiff + hadv + vdiff;
  end
end

flag = 0;
new_data = [];

return

% ===========================================================================

function [qd, flag, new_data] = quadfn(t, u, data)
% Quadrature integrand function

mx = data.mx;
my = data.my;
dx = data.dx;
dy = data.dy;
xmax = data.xmax;
ymax = data.ymax;

qd1 = 0.0;
for j = 1:my
  for i = 1:mx
    uij = u(j+(i-1)*my);
    if j == 1 | j == mx
      del_y = dy/2;
    else
      del_y = dy;
    end
    if i == 1 | i == mx
      del_x = dx/2;
    else
      del_x = dx;
    end
    qd1 = qd1 + uij * del_x*del_y;
  end
end

qd1 = qd1 / (xmax*ymax);

qd(1) = qd1;

flag = 0;
new_data = [];

return

% ===========================================================================

function [J, flag, new_data] = bjacfn(t, y, fy, data)
% Band Jacobian function

mx = data.mx;
my = data.my;
hordc = data.hdcoef;
horac = data.hacoef;
verdc = data.vdcoef;

mu = my;
ml = my;
mband = mu + 1 + ml;

for i = 1:mx
  for j = 1:my
     k = j + (i-1)*my;
     J(mu+1,k) = -2.0 * (verdc + hordc);
     if  i ~= 1
       J(1,k) = hordc + horac;
     end
     if i ~= mx
       J(mband,k) = hordc - horac;
     end
     if j ~= 1
       J(mu,k) = verdc;
     end
     if j ~= my
       J(mu+2,k) = verdc;
     end
  end
end

flag = 0;
new_data = [];
