function mkinFerTron_dns
% mkinFerTron_dns - Ferraris-Tronconi test problem
%
% Source: "Handbook of Test Problems in Local and Global Optimization",
%             C.A. Floudas, P.M. Pardalos et al.
%             Kluwer Academic Publishers, 1999.
% Test problem 4 from Section 14.1, Chapter 14: Ferraris and Tronconi
% 
% This problem involves a blend of trigonometric and exponential terms.
%    0.5 sin(x1 x2) - 0.25 x2/pi - 0.5 x1 = 0
%    (1-0.25/pi) ( exp(2 x1)-e ) + e x2 / pi - 2 e x1 = 0
% such that
%    0.25 <= x1 <=1.0
%    1.5 <= x2 <= 2 pi
% 
% The treatment of the bound constraints on x1 and x2 is done using
% the additional variables
%    l1 = x1 - x1_min >= 0
%    L1 = x1 - x1_max <= 0
%    l2 = x2 - x2_min >= 0
%    L2 = x2 - x2_max >= 0
% 
% and using the constraint feature in KINSOL to impose
%    l1 >= 0    l2 >= 0
%    L1 <= 0    L2 <= 0
% 
% The Ferraris-Tronconi test problem has two known solutions.
% The nonlinear system is solved by KINSOL using different 
% combinations of globalization and Jacobian update strategies 
% and with different initial guesses (leading to one or the other
% of the known solutions).
%
% Constraints are imposed to make all components of the solution
% positive.

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


% Initializations
% ---------------

% User data
lb = [0.25 ; 1.5];
ub = [1.0 ; 2*pi];
data.lb = lb;
data.ub = ub;

% Number of problem variables
nvar = 2;

% Number of equations = number of problem vars. + number of bound variables
neq = nvar + 2*nvar;

% Function tolerance
ftol = 1.0e-5;

% Step tolerance
stol = 1.0e-5;

% Constraints
constraints = [ 0 0 1 -1 1 -1];

% Modified/exact Newton
%   msbset = 0   -> modified Newton
%   msbset = 1   -> exact Newton
msbset = 0;

% Initialize solver
options = KINSetOptions('UserData', data, ...
                        'FuncNormTol', ftol, ...
                        'ScaledStepTol', stol, ...
                        'Constraints', constraints, ...
                        'MaxNumSetups', msbset, ...
                        'LinearSolver', 'Dense');
KINInit(@sysfn, neq, options);

% Initial guess
% -------------
%
% There are two known solutions for this problem
% 
% the following initial guess should take us to (0.29945; 2.83693)
x1 = lb(1);
x2 = lb(2);
u1 = [ x1 ; x2 ; x1-lb(1) ; x1-ub(1) ; x2-lb(2) ; x2-ub(2) ];
% while this one should take us to (0.5; 3.1415926)
x1 = 0.5*(lb(1)+ub(1));
x2 = 0.5*(lb(2)+ub(2));
u2 = [ x1 ; x2 ; x1-lb(1) ; x1-ub(1) ; x2-lb(2) ; x2-ub(2) ];

% No Y and F scaling
yscale = ones(neq,1);
fscale = ones(neq,1);

% No globalization
strategy = 'None';


fprintf('\nFerraris and Tronconi test problem\n');
fprintf('Tolerance parameters:\n');
fprintf('  fnormtol  = %10.6g\n  scsteptol = %10.6g\n', ftol, stol);
if msbset == 1
  fprintf('Exact Newton');
else
  fprintf('Modified Newton');
end
if strcmp(strategy,'None')
  fprintf('\n');
else
  fprintf(' with line search\n');
end

% Solve problem starting from the 1st initial guess
% -------------------------------------------------

fprintf('\n------------------------------------------\n');
fprintf('\nInitial guess on lower bounds\n');
fprintf('  [x1,x2] = %8.6g  %8.6g', u1(1), u1(2));

[status, u1] = KINSol(u1, strategy, yscale, fscale);
stats = KINGetStats;

fprintf('\nsolution\n');
fprintf('  [x1,x2] = %8.6g  %8.6g', u1(1), u1(2));
fprintf('\nSolver statistics:\n');
fprintf('  nni = %5d    nfe  = %5d \n', stats.nni, stats.nfe);
fprintf('  nje = %5d    nfeD = %5d \n', stats.LSInfo.njeD, stats.LSInfo.nfeD);

% Solve problem starting from the 2nd initial guess
% -------------------------------------------------

fprintf('\n------------------------------------------\n');
fprintf('\nInitial guess in middle of feasible region\n');
fprintf('  [x1,x2] = %8.6g  %8.6g', u2(1), u2(2));

[status, u2] = KINSol(u2, strategy, yscale, fscale);
stats = KINGetStats;

fprintf('\nsolution\n');
fprintf('  [x1,x2] = %8.6g  %8.6g', u2(1), u2(2));
fprintf('\nSolver statistics:\n');
fprintf('  nni = %5d    nfe  = %5d \n', stats.nni, stats.nfe);
fprintf('  nje = %5d    nfeD = %5d \n', stats.LSInfo.njeD, stats.LSInfo.nfeD);

% Free memory 
% --------------------------------------

KINFree;

return


% System function 
% ---------------

function [fu, flag, new_data] = sysfn(u, data)

lb = data.lb;
ub = data.ub;

x1 = u(1);
x2 = u(2);
l1 = u(3);
L1 = u(4);
l2 = u(5);
L2 = u(6);

e = exp(1);

fu(1) = 0.5 * sin(x1*x2) - 0.25 * x2 / pi - 0.5 * x1;
fu(2) = (1.0 - 0.25/pi)*(exp(2.0*x1)-e) + e*x2/pi - 2.0*e*x1;
fu(3) = l1 - x1 + lb(1);
fu(4) = L1 - x1 + ub(1);
fu(5) = l2 - x2 + lb(2);
fu(6) = L2 - x2 + ub(2);


flag = 0;
new_data = [];
return
