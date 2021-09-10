function midasBruss_ASA_dns
%midasBruss_ASA_dns - ASA example for the Brusselator problem
%  This example solves the forward and adjoint problems for
%  the 2D Brusselator example.
%
%  See also: midasBruss_dns

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

%--------------
% Problem data
%--------------

eps = 2.0e-3; % diffusion param
A = 1.0;
B = 3.4;

% Spatial length 0 <= x,y, <= L
L = 1.0;

% grid size
mx = 20;
my = 20;
dx = L/mx;
dy = L/my;

% coefficients in central FD
hdif = eps/dx^2;
vdif = eps/dy^2;

% problem dimension
nx = mx+1;
ny = my+1;
n = 2*nx*ny;

x = linspace(0,L,nx);
y = linspace(0,L,ny);

% Load user data structure
data.eps = eps;
data.A = A;
data.B = B;
data.L = L;
data.dx = dx;
data.dy = dy;
data.nx = nx;
data.ny = ny;
data.x = x;
data.y = y;
data.hdif = hdif;
data.vdif = vdif;

% ---------------------
% Initialize integrator
% ---------------------

% Integration limits

t0 = 0.0;
tf = 1.0;

% Initial conditions

[u, v] = BRUSic(data);
Y = UV2Y(u, v, data);
Yp = zeros(n,1);

% Specify algebraic variables

u_id = ones(ny,nx);
u_id(1,:) = 0;
u_id(ny,:) = 0;
u_id(:,1) = 0;
u_id(:,nx) = 0;

v_id = ones(ny,nx);
v_id(1,:) = 0;
v_id(ny,:) = 0;
v_id(:,1) = 0;
v_id(:,nx) = 0;

id = UV2Y(u_id, v_id, data);

% Optional inputs

options = IDASetOptions('UserData',data,...
                        'RelTol',1.e-5,...
                        'AbsTol',1.e-5,...
                        'VariableTypes',id,...
                        'suppressAlgVars','on',...
                        'MaxNumSteps', 1000,...
                        'LinearSolver','Dense');

% Initialize forward problem

IDAInit(@BRUSres,t0,Y,Yp,options);

% Compute consistent I.C.

[status, Y, Yp] = IDACalcIC(tf, 'FindAlgebraic');

% --------------
% Initialize ASA
% --------------

IDAAdjInit(150, 'Hermite');

% ---------------
% Integrate to tf
% ---------------

[status, t, Y] = IDASolve(tf, 'Normal');

[u,v] = Y2UV(Y, data);

stats_fwd = IDAGetStats;

plotSol(t,Y,data);

% ---------------------------
% Initialize backward problem
% ---------------------------

% Specify algebraic variables

l_id = ones(ny,nx);
l_id(1,:) = 0;
l_id(ny,:) = 0;
l_id(:,1) = 0;
l_id(:,nx) = 0;

m_id = ones(ny,nx);
m_id(1,:) = 0;
m_id(ny,:) = 0;
m_id(:,1) = 0;
m_id(:,nx) = 0;

idB = UV2Y(l_id, m_id, data);

% Final conditions

l = ones(ny,nx);
m = zeros(ny,nx);
YB = UV2Y(l, m, data);

lp = -2.0 * u .* v .* l + (B+1) * l;
mp = - l .* (u.^2);
YBp = UV2Y(lp, mp, data);

% Optional inputs

optionsB = IDASetOptions('UserData',data,...
                         'RelTol',1.e-5,...
                         'AbsTol',1.e-5,...
                         'VariableTypes',id,...
                         'suppressAlgVars','on',...
                         'LinearSolver','Dense');

% Initialize backward problem

idxB = IDAInitB(@BRUSresB,tf,YB,YBp,optionsB);

% --------------------------
% Backward integration to t0
% --------------------------

[status, t, YB] = IDASolveB(t0,'Normal');

plotSol(t,YB,data);

% -----------
% Free memory
% -----------

IDAFree;

return






% ====================================================================================
% Initial conditions
% ====================================================================================

function [u0, v0] = BRUSic(data)

dx = data.dx;
dy = data.dy;
nx = data.nx;
ny = data.ny;
L  = data.L;
x  = data.x;
y  = data.y;

n = 2*nx*ny;

[x2D , y2D] = meshgrid(x,y);

u0 = 1.0 - 0.5 * cos(pi*y2D/L);
u0(1,:) = u0(2,:);
u0(ny,:) = u0(ny-1,:);
u0(:,1) = u0(:,2);
u0(:,nx) = u0(:,nx-1);

v0 = 3.5 - 2.5*cos(pi*x2D/L);
v0(1,:) = v0(2,:);
v0(ny,:) = v0(ny-1,:);
v0(:,1) = v0(:,2);
v0(:,nx) = v0(:,nx-1);

return





% ====================================================================================
% Residual function
% ====================================================================================

function [res, flag, new_data] = BRUSres(t,Y,Yp,data)

nx = data.nx;
ny = data.ny;

A = data.A;
B = data.B;

hdif = data.hdif;
vdif = data.vdif;

% Convert Y and Yp to (u,v) and (up, vp)

[u,v] = Y2UV(Y,data);
[up,vp] = Y2UV(Yp,data);

% 2D residuals

ru = zeros(ny,nx);
rv = zeros(ny,nx);

% Inside the domain

for iy = 2:ny-1

  for ix = 2:nx-1
    
    uu = u(iy,ix);
    vv = v(iy,ix);
    
    ru(iy,ix) = up(iy,ix) ...
        - hdif * ( u(iy,ix+1) - 2*uu + u(iy,ix-1) ) ...
        - vdif * ( u(iy+1,ix) - 2*uu + u(iy-1,ix) ) ...
        - A + (B+1)*uu - uu^2 * vv;

    rv(iy,ix) = vp(iy,ix) ...
        - hdif * ( v(iy,ix+1) - 2*vv + v(iy,ix-1) ) ...
        - vdif * ( v(iy+1,ix) - 2*vv + v(iy-1,ix) ) ...
        - B*uu + uu^2 * vv;
    
  end
  
end
    
% Boundary conditions

ru(1,:) = u(1,:) - u(2,:);
ru(ny,:) = u(ny,:) - u(ny-1,:);
ru(:,1) = u(:,1) - u(:,2);
ru(:,nx) = u(:,nx) - u(:,nx-1);

rv(1,:) = v(1,:) - v(2,:);
rv(ny,:) = v(ny,:) - v(ny-1,:);
rv(:,1) = v(:,1) - v(:,2);
rv(:,nx) = v(:,nx) - v(:,nx-1);

% Convert (ru,rv) to res

res = UV2Y(ru,rv,data);

% Return flag and pb. data

flag = 0;
new_data = [];




% ====================================================================================
% Backward residual function
% ====================================================================================

function [resB, flag, new_data] = BRUSresB(t, Y, Yp, YB, YBp, data)

nx = data.nx;
ny = data.ny;

A = data.A;
B = data.B;

hdif = data.hdif;
vdif = data.vdif;

% Convert Y to (u,v)

[u,v] = Y2UV(Y,data);

% Convert YB and YBp to (l,m) and (lp,mp)

[l,m] = Y2UV(YB,data);
[lp,mp] = Y2UV(YBp,data);

% 2D residuals

rl = zeros(ny,nx);
rm = zeros(ny,nx);

% Inside the domain

for iy = 2:ny-1

  for ix = 2:nx-1

    uu = u(iy,ix);
    vv = v(iy,ix);
    
    ll = l(iy,ix);
    mm = m(iy,ix);
   
    rl(iy,ix) = lp(iy,ix) ...
        + hdif * ( l(iy,ix+1) - 2*ll + l(iy,ix-1) ) ...
        + vdif * ( l(iy+1,ix) - 2*ll + l(iy-1,ix) ) ...
        + 2*uu*vv*ll - (B+1)*ll + B*mm - 2*uu*vv*mm;

    rm(iy,ix) = mp(iy,ix) ...
        + hdif * ( m(iy,ix+1) - 2*mm + m(iy,ix-1) ) ...
        + vdif * ( m(iy+1,ix) - 2*mm + m(iy-1,ix) ) ...
        + ll * uu^2 - mm * uu^2;
    
  end
  
end

% Boundary conditions

rl(1,:) = l(1,:) - l(2,:);
rl(ny,:) = l(ny,:) - l(ny-1,:);
rl(:,1) = l(:,1) - l(:,2);
rl(:,nx) = l(:,nx) - l(:,nx-1);

rm(1,:) = m(1,:) - m(2,:);
rm(ny,:) = m(ny,:) - m(ny-1,:);
rm(:,1) = m(:,1) - m(:,2);
rm(:,nx) = m(:,nx) - m(:,nx-1);

% Convert (rl,rm) to resB

resB = UV2Y(rl,rm,data);

% Return flag and pb. data

flag = 0;
new_data = [];





% ====================================================================================
% 1D <-> 2D conversion functions
% ====================================================================================

function y = UV2Y(u, v, data)

nx = data.nx;
ny = data.ny;

u1 = reshape(u, 1, nx*ny);
v1 = reshape(v, 1, nx*ny);

y = reshape([u1;v1], 2*nx*ny,1);

return


function [u,v] = Y2UV(y, data)

nx = data.nx;
ny = data.ny;

y2 = reshape(y, 2, nx*ny);

u = reshape(y2(1,:), ny, nx);
v = reshape(y2(2,:), ny, nx);

return




% ====================================================================================
% Plot (u,v)
% ====================================================================================

function plotSol(t,Y,data)

x = data.x;
y = data.y;

[u,v] = Y2UV(Y, data);

figure;
set(gcf,'position',[600 600 650 300])

subplot(1,2,1)
surfc(x,y,u);
shading interp
%view(0,90)
view(-15,35)
axis tight
box on
grid off
xlabel('x');
ylabel('y');
title(sprintf('u(x,y,%g)',t))
%colorbar('horiz');

subplot(1,2,2)
surfc(x,y,v);
shading interp
%view(0,90)
view(-15,35)
axis tight
box on
grid off
xlabel('x');
ylabel('y');
title(sprintf('v(x,y,%g)',t))
%colorbar('horiz');

return
