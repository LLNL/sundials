function [x,y,u0_2, u_2] = mcvsDiurnal_kry
%mcvsDiurnal_kry - CVODES example problem (serial, Spgmr)
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

%------------------------
% SET USER DATA STRUCTURE
%------------------------

ns = 2;

mx = 50;
my = 20;

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
  y(jy) = ymin + (jy - 1) * dy;
end
for jx = 1:mx
  x(jx) = xmin + (jx - 1) * dx;
end

for jy = 1:my
  cy = (0.1 * (y(jy) - ymid))^2;
  cy = 1.0 - cy + 0.5 * cy^2;
  for jx = 1:mx
    cx = (0.1 * (x(jx) - xmid))^2;
    cx = 1.0 - cx + 0.5 * cx^2;
    u0(1,jx,jy) = c1s * cx * cy;
    u0(2,jx,jy) = c2s * cx * cy;
  end
end

u0_2 = squeeze(u0(2,:,:));

u0 = reshape(u0,2*mx*my,1);


%------------------------
% SET CVODES OPTIONS
%------------------------

% Tolerances
rtol = 1.0e-5;
atol = 1.0e-3;

options = CVodeSetOptions('UserData',data,...
                          'RelTol',rtol, 'AbsTol',atol,...
                          'LinearSolver','GMRES',...
                          'PrecType','Left',...
                          'PrecSetupFn',@psetfn,...
                          'PrecSolveFn',@psolfn);

%mondata = struct;
%options = CVodeSetOptions(options,'MonitorFn',@CVodeMonitor,'MonitorData',mondata);

CVodeInit(@rhsfn, 'BDF', 'Newton', t0, u0, options);

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

u_2 = squeeze(u(2,:,:));

si = CVodeGetStats

%------------------------
% FREE MEMORY
%------------------------

CVodeFree;

%------------------------
% PLOT RESULTS
%------------------------

figure;

hsurf0 = surf(y,x,u0_2);

set(hsurf0,'Edgecolor','flat','EdgeAlpha',0.6);
set(hsurf0,'FaceAlpha',0);

hold on

hsurf = surf(y,x,u_2);
set(hsurf,'FaceColor','interp', 'EdgeColor','none');

title('Initial and final values for species 2');
xlabel('y');
ylabel('x');
colorbar
box on

return

% ===========================================================================

function [ud, flag, new_data] = rhsfn(t, u, data)
% Right-hand side function

mx = data.mx;
xmin = data.xmin;
dx = data.dx;

my = data.my;
ymin = data.ymin;
dy = data.dy;

om = data.om;
q1 = data.q1;
q2 = data.q2;
c3 = data.c3;
a3 = data.a3;
a4 = data.a4;
hdco = data.hdco;
haco = data.haco;
vdco = data.vdco;

s = sin(om*t);
if s > 0.0
  q3 = exp(-a3/s);
  q4 = exp(-a4/s);
else
  q3 = 0.0;
  q4 = 0.0;
end

u = reshape(u, 2,mx*my);

for jy = 1:my
  ydn = ymin + (jy - 1.5)*dy;
  yup = ydn + dy;
  cydn = vdco * exp(0.2*ydn);
  cyup = vdco * exp(0.2*yup);
  i = (jy-1)*mx;
  idn = -mx;
  if jy == 1
    idn = mx;
  end
  iup = mx;
  if jy == my
    iup = -mx;
  end
  for jx = 1:mx
    ii = i + jx;
    c1 = u(1,ii);
    c2 = u(2,ii);
    % kinetic rate terms
    qq1 = q1 * c1 * c3;
    qq2 = q2 * c1 * c2;
    qq3 = q3 * c3;
    qq4 = q4 *c2;
    rkin1 = -qq1 - qq2 + 2.0*qq3 + qq4;
    rkin2 = qq1 - qq2 - qq4;
    % vertical diffusion
    c1dn = u(1,ii+idn);
    c2dn = u(2,ii+idn);
    c1up = u(1,ii+iup);
    c2up = u(2,ii+iup);
    vertd1 = cyup*(c1up - c1) - cydn*(c1 - c1dn);
    vertd2 = cyup*(c2up - c2) - cydn*(c2 - c2dn);
    % horizontal diffusion and advection
    ileft = -1;
    if jx == 1
      ileft = 1;
    end
    iright = 1;
    if jx == mx
      iright = -1;
    end
    c1lt = u(1,ii+ileft);
    c2lt = u(2,ii+ileft);
    c1rt = u(1,ii+iright);
    c2rt = u(2,ii+iright);
    hord1 = hdco * (c1rt-2.0*c1+c1lt);
    hord2 = hdco * (c2rt-2.0*c2+c2lt);
    horad1 = haco * (c1rt-c1lt);
    horad2 = haco * (c2rt-c2lt);
    % load into ud
    ud(1,ii) = vertd1 + hord1 + horad1 + rkin1; 
    ud(2,ii) = vertd2 + hord2 + horad2 + rkin2;
  end
  
end

ud = reshape(ud,2*mx*my,1);

flag = 0;
new_data = [];

new_data = data;

return

% ===========================================================================

function [jcur, flag, data] = psetfn(t,u,fu,jok,gm,data)
% Preconditioner setup function

persistent Jbd

mx = data.mx;
my = data.my;

if jok

  % Copy Jbd to P
  
  P = Jbd;
  jcur = false;
  
else

  % Generate Jbd from scratch and copy to P

  xmin = data.xmin;
  dx = data.dx;

  ymin = data.ymin;
  dy = data.dy;

  om = data.om;
  q1 = data.q1;
  q2 = data.q2;
  c3 = data.c3;
  a3 = data.a3;
  a4 = data.a4;
  hdco = data.hdco;
  haco = data.haco;
  vdco = data.vdco;
  
  s = sin(om*t);
  if s > 0.0
    q4 = exp(-a4/s);
  else
    q4 = 0.0;
  end

  u = reshape(u,2,mx*my);
  
  for jy = 1:my
    ydn = ymin + (jy - 1.5)*dy;
    yup = ydn + dy;
    cydn = vdco * exp(0.2*ydn);
    cyup = vdco * exp(0.2*yup);
    diag = -(cydn + cyup + 2.0*hdco);
    i = (jy-1)*mx;
    for jx = 1:mx
      ii = i + jx;
      c1 = u(1,ii);
      c2 = u(2,ii);
      Jbd(1,1,ii) = (-q1*c3 - q2*c2) + diag;
      Jbd(1,2,ii) = -q2*c1 + q4;
      Jbd(2,1,ii) = q1*c3 - q2*c2;
      Jbd(2,2,ii) = (-q2*c1 - q4) + diag;
    end
  end
  
  P = Jbd;
  jcur = true;

end

% Scale by -gamma and add identity
P = - gm*P;
for i = 1:mx*my
  P(:,:,i) = eye(2) + P(:,:,i);
end

flag = 0;
data.P = P;


return

% ===========================================================================

function [z, flag, new_data] = psolfn(t,y,fy,r,data)
% Preconditioner solve function

P = data.P;
 
mx = data.mx;
my = data.my;

r = reshape(r,2,mx*my);

for i = 1:mx*my
  z(:,i) = P(:,:,i)^(-1)*r(:,i);
end

z = reshape(z,2*mx*my,1);

flag = 0;
new_data = [];

return