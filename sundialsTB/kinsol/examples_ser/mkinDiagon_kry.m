function mkinDiagon_kry
%mkinDiagon_kry - KINSOL example problem (serial, GMRES)
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
%   See also: kindiag_sys kindag_pset kindiag_psol

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
% $Revision$Date: 2007/10/26 16:30:49 $

neq = 128;

strategy = 'None';

fnormtol = 1.0e-5;
scsteptol = 1.0e-4;
maxl   = 10;
maxrs  = 2;
msbset = 5;

data.P = [];

options = KINSetOptions('UserData', data,...
                        'Verbose',true,...
                        'FuncNormTol', fnormtol,...
                        'ScaledStepTol',scsteptol,...
                        'LinearSolver','GMRES',....
                        'KrylovMaxDim', maxl,...
                        'MaxNumRestarts', maxrs,...
                        'MaxNumSetups', msbset,...
                        'PrecSetupFn',@psetfn,...
                        'PrecSolveFn',@psolfn);
KINInit(@sysfn, neq, options);

y0 = 2.0*[1:neq]';
scale = ones(neq,1);

[status, y] = KINSol(y0, strategy, scale, scale);

if status < 0
  fprintf('KINSOL failed. status = %d\n',status);
else
  for i = 1:4:neq
    fprintf('%4d   |  %6.2f  %6.2f  %6.2f  %6.2f\n',...
            i, y(i), y(i+1), y(i+2), y(i+3));
  end
end

stats = KINGetStats;
ls_stats = stats.LSInfo;

stats
ls_stats

KINFree;


% ============================================================

function [fy, flag, new_data] = sysfn(y, data)

neq = length(y);
for i = 1:neq
  fy(i) = y(i)^2 - i^2;
end

new_data = [];  % data was not modified
flag = 0;       % success

% ============================================================

function [flag, new_data] = psetfn(y,yscale,fy,fscale,data)

neq = length(y);

for i = 1:neq
  P(i) = 0.5 / (y(i)+5.0);
end

new_data.P = P;  % updated P in data structure
flag = 0;        % success

% ============================================================

function [x, flag, new_data] = psolfn(y,yscale,fy,fscale,v,data)

P = data.P;

neq = length(y);

for i=1:neq
  x(i) = v(i) * P(i);
end

new_data = []; % data was not modified
flag = 0;      % success
