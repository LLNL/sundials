function [] = mkinDiagon_kry_p(comm)
%mkinDiagon_kry_p - KINSOL example problem (parallel, GMRES)
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
%   See also: mpirun

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

% User data structure

data.mype = mype;       % MPI id
data.nlocal = nlocal;   % local problem size
data.P = [];            % workspace for preconditioner


options = KINSetOptions('UserData', data,...
                        'FuncNormTol', fnormtol,...
                        'ScaledStepTol',scsteptol,...
                        'LinearSolver','GMRES',....
                        'KrylovMaxDim', maxl,...
                        'MaxNumRestarts', maxlrst,...
                        'MaxNumSetups', msbset,...
                        'PrecSetupFn',@psetfn,...
                        'PrecSolveFn',@psolfn);

if mype==0
  options = KINSetOptions(options,'Verbose',true);
end

KINInit(@sysfn, nlocal, options);

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

% =======================================================

function [fy, flag, new_data] = sysfn(y, data)

nlocal = data.nlocal;
mype = data.mype;
baseadd = mype * nlocal;

for i = 1:nlocal
  fy(i) = y(i)^2 - (i+baseadd)^2;
end

new_data = []; % data was not modified
flag = 0;      % success

% =======================================================

function [flag, new_data] = psetfn(y,yscale,fy,fscale,data)

nlocal = data.nlocal;

for i = 1:nlocal
  P(i) = 0.5 / (y(i)+5.0);
end

new_data = data;
new_data.P = P;

flag = 0;

% =======================================================

function [x, flag, new_data] = psolfn(y,yscale,fy,fscale,v,data)

nlocal = data.nlocal;
P = data.P;

for i=1:nlocal
  x(i) = v(i) * P(i);
end

new_data = [];
flag = 0;
