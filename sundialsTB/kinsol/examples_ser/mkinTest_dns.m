function mkinTest_dns
%mkinTest_dns - KINSOL example problem (serial, dense)
%   Simple test problem for the Dense linear solver in KINSOL
%   This example solves the system
%       y(1)^2 + y(2)^2 = 1
%       y(2) = y(1)^2
%

% Radu Serban <radu@llnl.gov>
% LLNS Start Copyright
% Copyright (c) 2013, Lawrence Livermore National Security
% This work was performed under the auspices of the U.S. Department 
% of Energy by Lawrence Livermore National Laboratory in part under 
% Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
% Produced at the Lawrence Livermore National Laboratory.
% All rights reserved.
% For details, see the LICENSE file.
% LLNS End Copyright
% $Revision: 1.2 $Date: 2007/10/26 16:30:49 $

% Initialize problem
neq = 2;
fnormtol  = 1.0e-5;
scsteptol = 1.0e-4;
msbset = 1; % force exact Newton
options = KINSetOptions('FuncNormTol', fnormtol,...
                        'ScaledStepTol',scsteptol,...
                        'LinearSolver','Dense',....
                        'MaxNumSetups', msbset);
KINInit(@sysfn, neq, options);

% Solve problem
y0 = ones(neq,1);
scale = ones(neq,1);
strategy = 'LineSearch';
[status, y] = KINSol(y0, strategy, scale, scale);

% Evaluate system function at solution
[fy, flag] = sysfn(y);

% Print results
fprintf('Solution: %10.4e  %10.4e\n', y(1), y(2));
fprintf('Residual: %10.4e  %10.4e\n', fy(1), fy(2));

slv_stats = KINGetStats;
ls_stats = slv_stats.LSInfo;
slv_stats
ls_stats


% Free memory
KINFree;

return

% ===================================================================

function [fy, flag] = sysfn(y)

fy(1) = y(1)^2 + y(2)^2 - 1.0;
fy(2) = y(2) - y(1)^2;

flag = 0;

return


