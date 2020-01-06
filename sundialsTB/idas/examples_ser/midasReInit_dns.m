function [] = midasReInit_dns()
%midasReInit_dns - Illustration of the IDAS reinitialization function

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


fprintf('Example for integrating over a discontinuity in states\n');
fprintf('using the IDAS re-initialization function\n\n');  
fprintf('Integrate over t = [ 0 1.0 ] the DAE:');
fprintf('  y1'' + y1 - y2 = 0\n');
fprintf('         y1 + y2 = 0\n');
fprintf('with initial conditions:\n');
fprintf('  y1(0) =  1.0\n');
fprintf('  y2(0) = -1.0\n');
fprintf('until y2(t*) = -0.5. At t*, perturb:\n');
fprintf('  y1(t*) <- y1(t*) - 0.25\n');
fprintf('  y2(t*) <- y2(t*) + 0.25\n');
fprintf('and continue the integration to t = 1.0\n\n');

t0 = 0.0;
tout = 1.0;
y0 = [1.0;-1.0];
yp0 = [-2.0;0.0];

% Set optional inputs
options = IDASetOptions('RelTol',1.e-4,...
                        'AbsTol',1.e-5,...
                        'LinearSolver','Dense');
options = IDASetOptions(options,'RootsFn',@my_rootfct, 'NumRoots',1);

% Initialize solver
IDAInit(@my_resfct,t0,y0,yp0,options);

% Initialize arrays
tt = [];
yy1 = [];
yy2 = [];

% Integrate DAE until root is found
t = t0;
while t<tout
  [status, t, y] = IDASolve(tout,'OneStep');
  tt = [tt;t];
  yy1 = [yy1;y(1)];
  yy2 = [yy2;y(2)];
  if status == 2
    break;
  end
end


fprintf('');

% Get yp at current time
yp = IDAGet('DerivSolution',t,0);

% Add discontinuity in solution 
% (must be consistent with the algebraic constraint)
t0 = t;
y0 = [y(1)-0.25; y(2)+0.25];
yp0 = [yp(1)+0.25; 0.0];

% Reinitialize solver
IDAReInit(t0,y0,yp0,options);

% Integrate to final time
t = t0;
while t<tout
  [status, t, y] = IDASolve(tout,'OneStep');
  tt = [tt;t];
  yy1 = [yy1;y(1)];
  yy2 = [yy2;y(2)];
end

% Free memory
IDAFree;

% Plot solution
hp1 = plot(tt,yy1);
hold on
hp2 = plot(tt,yy2,'r');
set(gca,'XLim',[tt(1) tt(end)]);
plot([t0 t0],get(gca,'YLim'),'k:');

legend([hp1,hp2],'y_1','y_2');

%
% ========================================================
%

function [rr, flag] = my_resfct(t, y, yp)

rr(1) = yp(1) + y(1) - y(2);
rr(2) = y(1) + y(2);

flag = 0;

%
% ========================================================
%

function [g, flag] = my_rootfct(t, y, yp)

g(1) = y(2) + 0.5;

flag = 0;
