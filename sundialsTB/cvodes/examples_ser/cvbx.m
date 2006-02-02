%CVBX - CVODES example problem (serial, band)
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
%
%   See also: cvbx_f, cvbx_q, cvbx_J

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/01/06 18:59:48 $

xmax = 2.0;
ymax = 1.0;
mx = 10;
my = 5;

rtol = 0.0;
atol = 1.0e-5;

t0 = 0.0;
dtout = 0.1;
nout = 10;

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

% initial conditions for states
t = t0;
u = zeros(mx*my,1);
for j = 1:my
  y = j * data.dy;
  for i = 1:mx
    x = i * data.dx;
    u(j + (i-1)*my) = x*(xmax-x)*y*(ymax-y)*exp(5.0*x*y);
  end
end

% Initial condition for quadrature variable
q = 0.0;

options = CVodeSetOptions('RelTol',rtol, 'AbsTol',atol);
options = CVodeSetOptions(options,...
                          'LinearSolver','Band',...
                          'JacobianFn',@cvbx_J,...
                          'UpperBwidth',my,...
                          'LowerBwidth',my);
options = CVodeSetOptions(options,...
                          'Quadratures', 'on',...
                          'QuadRhsFn',@cvbx_q,...
                          'QuadInitCond', q);




%options = CVodeSetOptions(options,'MonitorFn','CVodeMonitor');

CVodeMalloc(@cvbx_f, t, u, options, data);


ff=figure;
hold on;
box

umax = norm(u,'inf');
uavg = cvbx_q(t,u,data);
fprintf('At t = %f   max.norm(u) = %e\n',t, umax);

for i = 1:nout

  t_old = t;
  uavg_old = uavg;

  tout = t + dtout;
  [status,t,u, q] = CVode(tout, 'Normal');
%  [status,t,u] = CVode(tout, 'Normal');

  if status ~= 0
    return
  end
  
  uavg = cvbx_q(t,u,data);
  umax = norm(u,'inf');
  fprintf('At t = %f   max.norm(u) = %e\n',tout, umax);

  figure(ff);
  plot([t_old t],[uavg_old uavg]);
%  plot([t0 tout], [q q]/(tout-t0), 'r-');
%  plot([tout tout], [0 q]/(tout-t0), 'r-');
  
end
  
si= CVodeGetStats

CVodeFree;