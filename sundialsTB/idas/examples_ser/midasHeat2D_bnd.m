function midasHeat2D_bnd
%midasHeat2D_bnd: 2D heat equation, serial, banded.
%
% This example solves a discretized 2D heat equation problem.
% This version uses the band solver IDABand, and IDACalcIC.
%
% The DAE system solved is a spatial discretization of the PDE
%          du/dt = d^2u/dx^2 + d^2u/dy^2
% on the unit square. The boundary condition is u = 0 on all edges.
% Initial conditions are given by u = 16 x (1 - x) y (1 - y).
% The PDE is treated with central differences on a uniform M x M
% grid. The values of u at the interior points satisfy ODEs, and
% equations u = 0 at the boundaries are appended, to form a DAE
% system of size N = M^2. Here M = 10.
%
% The system is solved with IDA using the banded linear system
% solver, half-bandwidths equal to M, and default
% difference-quotient Jacobian. For purposes of illustration,
% IDACalcIC is called to compute correct values at the boundary,
% given incorrect values as input initial guesses. The constraints
% u >= 0 are posed for all components. Output is taken at
% t = 0, .01, .02, .04, ..., 10.24. (Output at t = 0 is for
% IDACalcIC cost statistics only.)

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
% $Revision$Date: 2007/08/21 17:38:43 $

m = 20;
N = m^2;
data.m = m;
data.N = N;
data.dx = 1.0/(m-1);
data.c = 1.0/data.dx^2;

fp = figure;
set(gcf,'position',[250 175 560 900]);

[t0,yy0,yp0,id,cnstr] = ic(data); 

% Plot initial guess for IC

figure(fp);
subplot(2,1,1);
hold on
hs1 = surf(reshape(yy0,m,m));
shading interp
set(hs1,'FaceAlpha',0.35);
box on
view(-30,30)

options = IDASetOptions('UserData',data,...
                        'RelTol',0.0,...
                        'AbsTol',1.0e-3,...
                        'VariableTypes',id,...
                        'ConstraintTypes',cnstr,...
                        'LinearSolver','Band',...
                        'LowerBwidth',m,...
                        'UpperBwidth',m);

IDAInit(@resfun,t0,yy0,yp0,options);

tout = 0.01;
[status, yy0_mod, yp0_mod] = IDACalcIC(tout, 'FindAlgebraic');

% Plot corrected IC

figure(fp);
subplot(2,1,1);
hs1 = surf(reshape(yy0_mod,m,m));
set(hs1,'FaceColor','none');


% Plot solution

subplot(2,1,2);
hold on
hs1 = surf(reshape(yy0_mod,m,m));
shading interp
view(-30,30)
zlim_yy = get(gca,'ZLim');
box on

fprintf('t = %.4f    [Press any key]\n',t0);
pause;

nout = 5;
tout = 0.01;

for iout = 1:nout
  [status,t,yy] = IDASolve(tout,'Normal');
  tout = 2*tout;

  figure(fp);
  subplot(2,1,2);
  set(hs1,'FaceAlpha',0.15);
  hs1 = surf(reshape(yy,m,m));
  shading interp
  set(gca,'ZLim',zlim_yy);

  fprintf('t = %.4f    [Press any key]\n',t);
  pause;
  
end

IDAFree;


function [t,yy,yp,id,cnstr] = ic(data)

m = data.m;
N = data.N;
dx = data.dx;

id = ones(N,1);
cnstr = ones(N,1);
yy = zeros(N,1);
yp = zeros(N,1);

t = 0.0;

% Initialize yy on all grid points. */ 
for j=0:m-1
  yfact = dx * j;
  offset = m*j;
  for i=0:m-1
    xfact = dx * i;
    loc = offset + i + 1;
    yy(loc) = 16.0 * xfact * (1.0 - xfact) * yfact * (1.0 - yfact);
  end
end
  
% The residual gives the negative of ODE RHS values at 
% interior points.
yp = zeros(N,1);
[yp,flag,new_data] = resfun(t,yy,yp,data);
yp = -yp;

% Finally, set values of yy, yp, and id at boundary points.
for j=0:m-1
  offset = m*j;
  for i=0:m-1
    loc = offset + i + 1;
    if (j == 0 || j == m-1 || i == 0 || i == m-1 )
        yy(loc) = 0.1;
        yp(loc) = 0.0;
        id(loc) = 0.0;
    end
  end
end

% ====================================================================

function [rr,flag,new_data] = resfun(t,yy,yp,data)

m = data.m;
N = data.N;
dx = data.dx;
c = data.c;

% Initialize resval to uu, to take care of boundary equations.
rr = yy;
  
% Loop over interior points; set rr = yp - (central difference).
for j = 1:m-2
  offset = m*j;
  for i = 1:m-2
    loc = offset + i + 1;
    rr(loc) = yp(loc) - c * ...
	  (yy(loc-1) + yy(loc+1) + yy(loc-m) + yy(loc+m) - 4.0*yy(loc));
  end
end

flag = 0;
new_data = [];