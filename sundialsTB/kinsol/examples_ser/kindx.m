%KINDX - KINSOL example problem (serial, dense)
%   Simple test problem for the Dense linear solver in KINSOL
%   This example solves the system
%       y(1)^2 + y(2)^2 = 1
%       y(2) = y(1)^2
%
%   See also: kindx_sys

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$

neq = 2;

strategy = 'LineSearch';

fnormtol  = 1.0e-5;
scsteptol = 1.0e-4;
msbset = 1; % force exact Newton

options = KINSetOptions('FuncNormTol', fnormtol,...
                        'ScaledStepTol',scsteptol,...
                        'LinearSolver','Dense',....
                        'MaxNumSetups', msbset);

KINMalloc('kindx_sys', neq, options);

y0 = ones(neq,1);
scale = ones(neq,1);

[status, y] = KINSol(y0, strategy, scale, scale)

KINFree;