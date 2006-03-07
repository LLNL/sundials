%CVKX - CVODES example problem (serial, Spgmr)
%   An ODE system is generated from the following 2-species diurnal
%   kinetics advection-diffusion PDE system in 2 space dimensions:
%
%   dc(i)/dt = Kh*(d/dx)^2 c(i) + V*dc(i)/dx + (d/dy)(Kv(y)*dc(i)/dy)
%                   + Ri(c1,c2,t)      for i = 1,2,   where
%     R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 ,
%     R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 ,
%     Kv(y) = Kv0*exp(y/5) ,
%   Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t)
%   vary diurnally. The problem is posed on the square
%      0 <= x <= 20,    30 <= y <= 50   (all in km),
%   with homogeneous Neumann boundary conditions, and for time t in
%      0 <= t <= 86400 sec (1 day).
%   The PDE system is treated by central differences on a uniform
%   10 x 10 mesh, with simple polynomial initial profiles.
%   The problem is solved with CVODES, with the BDF/GMRES
%   method (i.e. using the CVSPGMR linear solver) and the
%   block-diagonal part of the Newton matrix as a left
%   preconditioner. A copy of the block-diagonal part of the
%   Jacobian is saved and conditionally reused within the Precond
%   routine.
%
%   See also: cvkx_f, cvkx_pset, cvkx_psol

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/01/06 18:59:49 $

%------------------------
% SET USER DATA STRUCTURE
%------------------------

ns = 2;

mx = 10;
my = 10;

xmin = 0.0;  xmax = 20.0; xmid = 10.0;
ymin = 30.0; ymax = 50.0; ymid = 40.0;
dx = (xmax-xmin)/(mx-1);
dy = (ymax-ymin)/(my-1);

kh = 4.0e-6;
vel = 0.001;
kv0 = 1.0e-8;
halfday = 4.32e4;

c1s = 1.0e6;
c2s = 1.0e12;


data.ns = ns;

%% Problem constants

data.mx = mx;
data.xmin = xmin;
data.xmax = xmax;
data.xmid = xmid;
data.dx = dx;

data.my = my;
data.ymin = ymin;
data.ymax = ymax;
data.ymid = ymid;
data.dy = dy;

data.q1 = 1.63e-16;
data.q2 = 4.66e-16;
data.c3 = 3.7e16;
data.a3 = 22.62;
data.a4 = 7.601;

data.om = pi/halfday;
data.hdco = kh/dx^2;
data.haco = vel/(2*dx);
data.vdco = kv0/dy^2;

%% Workspace

data.P = [];

%------------------------
% SET INITIAL PROFILE
%------------------------

t0 = 0.0;

for jy = 1:my
  y = ymin + (jy - 1) * dy;
  cy = (0.1 * (y - ymid))^2;
  cy = 1.0 - cy + 0.5 * cy^2;
  for jx = 1:mx
    x = xmin + (jx - 1) * dx;
    cx = (0.1 * (x - xmid))^2;
    cx = 1.0 - cx + 0.5 * cx^2;
    u0(1,jx,jy) = c1s * cx * cy;
    u0(2,jx,jy) = c2s * cx * cy;
  end
end

u0 = reshape(u0,2*mx*my,1);

%------------------------
% SET CVODES OPTIONS
%------------------------

% Tolerances
rtol = 1.0e-5;
atol = 1.0e-3;

options = CVodeSetOptions('RelTol',rtol, 'AbsTol',atol,...
                          'LinearSolver','GMRES',...
                          'PrecType','Left',...
                          'PrecSetupFn',@cvkx_pset,...
                          'PrecSolveFn',@cvkx_psol);

mondata = struct;
options = CVodeSetOptions(options,'MonitorFn',@CVodeMonitor,'MonitorData',mondata);

CVodeMalloc(@cvkx_f,t0,u0,options,data);

%------------------------
% SOLVE PROBLEM
%------------------------

twohr = 7200.0;
tout = twohr;
nout = 12;

for i = 1:nout
  [status,t,u] = CVode(tout,'Normal');
  si = CVodeGetStats;
  u = reshape(u,2,mx,my);
  fprintf('status = %d   t = %.2e   nst = %d   q = %d   h = %.2e\n',...
          status, t, si.nst, si.qlast, si.hlast);
  fprintf('c1 (bot.left/middle/top rt.) = %12.3e  %12.3e  %12.3e\n',...
          u(1,1,1), u(1,5,5), u(1,10,10));
  fprintf('c2 (bot.left/middle/top rt.) = %12.3e  %12.3e  %12.3e\n',...
          u(2,1,1), u(2,5,5), u(2,10,10));
  tout = tout + twohr;
end

si = CVodeGetStats

%------------------------
% FREE MEMORY
%------------------------

CVodeFree;


