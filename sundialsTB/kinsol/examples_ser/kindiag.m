%KINDIAG - KINSOL example problem (serial, GMRES)
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
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$

neq = 128;

strategy = 'None';

fnormtol = 1.0e-5;
scsteptol = 1.0e-4;
maxl   = 10;
maxrs  = 2;
msbset = 5;

options = KINSetOptions('FuncNormTol', fnormtol,...
                        'ScaledStepTol',scsteptol,...
                        'LinearSolver','GMRES',....
                        'KrylovMaxDim', maxl,...
                        'MaxNumRestarts', maxrs,...
                        'MaxNumSetups', msbset,...
                        'PrecSetupFn','kindiag_pset',...
                        'PrecSolveFn','kindiag_psol');
data.P = [];

KINMalloc('kindiag_sys', neq, options, data);

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

s = KINGetStats


KINFree;

