function [] = kindiagp(comm)
%KINDIAGP - KINSOL example problem (parallel, GMRES)
%   Simple diagonal test, using user-supplied preconditioner setup and 
%   solve routines.
%
%   This example does a basic test of the solver by solving the system:
%               f(y) = 0  for
%               f(y) = y(i)^2 - i^2
%
%   No scaling is done.
%   An approximate diagonal preconditioner is used.
%
%   See also: mpirun kindiagp_sys kindagp_pset kindiagp_psol

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/01/06 19:00:04 $

[status npes] = MPI_Comm_size(comm);
[status mype] = MPI_Comm_rank(comm);

% Problem size

nlocal = 32;
neq = npes * nlocal;

% Problem options

fnormtol = 1.0e-5;
scsteptol = 1.0e-4;
maxl = 10;
maxlrst = 2;
msbset = 5;

options = KINSetOptions('FuncNormTol', fnormtol,...
                        'ScaledStepTol',scsteptol,...
                        'LinearSolver','GMRES',....
                        'KrylovMaxDim', maxl,...
                        'MaxNumRestarts', maxlrst,...
                        'MaxNumSetups', msbset,...
                        'PrecSetupFn',@kindiagp_pset,...
                        'PrecSolveFn',@kindiagp_psol);

if mype==0
  options = KINSetOptions(options,'Verbose',true);
end

% User data structure

data.mype = mype;       % MPI id
data.nlocal = nlocal;   % local problem size
data.P = [];            % workspace for preconditioner

KINMalloc(@kindiagp_sys, nlocal, options, data);

% Initial guess and scale vector

y0 = 2.0 * ([1:nlocal] + mype*nlocal);
scale = ones(nlocal,1);

% Solve problem using Inexact Newton 

strategy = 'None';
[status, y] = KINSol(y0, strategy, scale, scale);

% Print solution

if status < 0
  fprintf('KINSOL failed. status = %d\n',status);
else
  for i = 1:4:nlocal
    fprintf('%4d   |  %6.2f  %6.2f  %6.2f  %6.2f\n',...
            i, y(i), y(i+1), y(i+2), y(i+3));
  end
end

fprintf('DONE\n');


% Free memory

KINFree;
