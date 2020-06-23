function mkinRoboKin_dns
% mkinRoboKin_dns - nonlinear system from robot kinematics.
%
% Source: "Handbook of Test Problems in Local and Global Optimization",
%             C.A. Floudas, P.M. Pardalos et al.
%             Kluwer Academic Publishers, 1999.
% Test problem 6 from Section 14.1, Chapter 14
% 
% The nonlinear system is solved by KINSOL using the DENSE linear solver.
%
% Constraints are imposed to make all components of the solution be
% within [-1,1]. This is achieved by introducing additional "bound variables",
%    l_i = x_i + 1   and u-i = 1 - x_i, i = 1,...,nvar
% and using the Constraints option to KINSOL to enforce l_i >=0 and u_i >= 0.
 
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
% $Revision$Date: 2007/10/26 16:30:49 $


fprintf('\nRobot Kinematics Example\n');
fprintf('8 variables; -1 <= x_i <= 1\n');
fprintf('KINSOL problem size: 8 + 2*8 = 24 \n\n');

% Number of problem variables
nvar = 8;

% Number of equations = number of problem vars. + number of bound variables
neq = nvar + 2*nvar;
 
% Function tolerance
ftol = 1.0e-5;

% Step tolerance
stol = 1.0e-5;

% Constraints (all bound variables non-negative)
constraints = [ zeros(nvar,1); ones(nvar,1) ; ones(nvar,1) ];

% Force exact Newton
msbset = 1;

% Initialize solver
options = KINSetOptions('FuncNormTol', ftol, ...
                        'ScaledStepTol', stol, ...
                        'Constraints', constraints, ...
                        'MaxNumSetups', msbset, ...
                        'LinearSolver', 'Dense', ...
                        'JacobianFn', @sysjac);
KINInit(@sysfn, neq, options);

% Initial guess
y = ones(neq,1);
y(1:nvar) = 1.0/sqrt(2);

fprintf('Initial guess:\n');
PrintOutput(y);

% Call KINSol to solve the problem
yscale = ones(neq,1);
fscale = ones(neq,1);
strategy = 'LineSearch';
[status, y] = KINSol(y, strategy, yscale, fscale);

fprintf('\nComputed solution:\n');
PrintOutput(y);

% Print final statistics and free memory 
stats = KINGetStats;
ls_stats = stats.LSInfo;
fprintf('\nSolver statistics:');
stats
ls_stats

KINFree;

return;

% System function 
% ---------------

function [fy, flag] = sysfn(y)

% Extract problem variables and bound variables

x1 = y(1); l1 = y( 9); u1 = y(17); 
x2 = y(2); l2 = y(10); u2 = y(18); 
x3 = y(3); l3 = y(11); u3 = y(19); 
x4 = y(4); l4 = y(12); u4 = y(20); 
x5 = y(5); l5 = y(13); u5 = y(21); 
x6 = y(6); l6 = y(14); u6 = y(22); 
x7 = y(7); l7 = y(15); u7 = y(23); 
x8 = y(8); l8 = y(16); u8 = y(24); 

% Nonlinear equations 

eq1 = - 0.1238*x1 + x7 - 0.001637*x2 - 0.9338*x4 + 0.004731*x1*x3 - 0.3578*x2*x3 - 0.3571;
eq2 = 0.2638*x1 - x7 - 0.07745*x2 - 0.6734*x4 + 0.2238*x1*x3 + 0.7623*x2*x3 - 0.6022;
eq3 = 0.3578*x1 + 0.004731*x2 + x6*x8;
eq4 = - 0.7623*x1 + 0.2238*x2 + 0.3461;
eq5 = x1*x1 + x2*x2 - 1;
eq6 = x3*x3 + x4*x4 - 1;
eq7 = x5*x5 + x6*x6 - 1;
eq8 = x7*x7 + x8*x8 - 1;

% Lower bounds ( l_i = 1 + x_i >= 0)

lb1 = l1 - 1.0 - x1;
lb2 = l2 - 1.0 - x2;
lb3 = l3 - 1.0 - x3;
lb4 = l4 - 1.0 - x4;
lb5 = l5 - 1.0 - x5;
lb6 = l6 - 1.0 - x6;
lb7 = l7 - 1.0 - x7;
lb8 = l8 - 1.0 - x8;

% Upper bounds ( u_i = 1 - x_i >= 0)

ub1 = u1 - 1.0 + x1;
ub2 = u2 - 1.0 + x2;
ub3 = u3 - 1.0 + x3;
ub4 = u4 - 1.0 + x4;
ub5 = u5 - 1.0 + x5;
ub6 = u6 - 1.0 + x6;
ub7 = u7 - 1.0 + x7;
ub8 = u8 - 1.0 + x8;

% Load residuals for the problem equations
% and the equations encoding the constraints

fy(1) = eq1; fy( 9) = lb1; fy(17) = ub1;
fy(2) = eq2; fy(10) = lb2; fy(18) = ub2;
fy(3) = eq3; fy(11) = lb3; fy(19) = ub3;
fy(4) = eq4; fy(12) = lb4; fy(20) = ub4;
fy(5) = eq5; fy(13) = lb5; fy(21) = ub5;
fy(6) = eq6; fy(14) = lb6; fy(22) = ub6;
fy(7) = eq7; fy(15) = lb7; fy(23) = ub7;
fy(8) = eq8; fy(16) = lb8; fy(24) = ub8;

flag = 0;

return;


% System Jacobian
% ---------------

function [J, flag] = sysjac(y, fy)

% Extract problem variables

x1 = y(1);
x2 = y(2);
x3 = y(3);
x4 = y(4);
x5 = y(5);
x6 = y(6);
x7 = y(7);
x8 = y(8);

% Nonlinear equations

%  - 0.1238*x1 + x7 - 0.001637*x2 - 0.9338*x4 + 0.004731*x1*x3 - 0.3578*x2*x3 - 0.3571 

J(1,1) = - 0.1238 + 0.004731*x3;
J(1,2) = - 0.001637 - 0.3578*x3;
J(1,3) = 0.004731*x1 - 0.3578*x2;
J(1,4) = - 0.9338;
J(1,7) = 1.0;

% 0.2638*x1 - x7 - 0.07745*x2 - 0.6734*x4 + 0.2238*x1*x3 + 0.7623*x2*x3 - 0.6022

J(2,1) = 0.2638 + 0.2238*x3;
J(2,2) = - 0.07745 + 0.7623*x3;
J(2,3) = 0.2238*x1 + 0.7623*x2;
J(2,4) = - 0.6734;
J(2,7) = -1.0;

% 0.3578*x1 + 0.004731*x2 + x6*x8

J(3,1) = 0.3578;
J(3,2) = 0.004731;
J(3,6) = x8;
J(3,8) = x6;

% - 0.7623*x1 + 0.2238*x2 + 0.3461

J(4,1) = - 0.7623;
J(4,2) = 0.2238;

% x1*x1 + x2*x2 - 1

J(5,1) = 2.0*x1;
J(5,2) = 2.0*x2;

% x3*x3 + x4*x4 - 1

J(6,3) = 2.0*x3;
J(6,4) = 2.0*x4;

% x5*x5 + x6*x6 - 1

J(7,5) = 2.0*x5;
J(7,6) = 2.0*x6;

% x7*x7 + x8*x8 - 1

J(8,7) = 2.0*x7;
J(8,8) = 2.0*x8;

% Lower bounds ( l_i = 1 + x_i >= 0)
%     l_i - 1.0 - x_i

for i=1:8
  J(8+i,i)   = -1.0;
  J(8+i,8+i) =  1.0;
end

% Upper bounds ( u_i = 1 - x_i >= 0)
%     u_i - 1.0 + x_i

for i=1:8
  J(16+i,i)    = 1.0;
  J(16+i,16+i) = 1.0;
end

flag = 0;

return



% Print solution
% --------------

function PrintOutput(y)

nvar = 8;

fprintf('     l=x+1          x         u=1-x\n');
fprintf('   ----------------------------------\n');

for i = 1:nvar
  fprintf(' %10.6g   %10.6g   %10.6g\n', y(i+nvar), y(i), y(i+2*nvar));
end

return

