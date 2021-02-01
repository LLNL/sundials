function [status, y] = KINSol(y0, strategy, yscale, fscale)
%KINSol solves the nonlinear problem.
%
%   Usage:  [STATUS, Y] = KINSol(Y0, STRATEGY, YSCALE, FSCALE)
%
%   KINSol manages the computational process of computing an approximate 
%   solution of the nonlinear system. If the initial guess (initial value 
%   assigned to vector Y0) doesn't violate any user-defined constraints, 
%   then KINSol attempts to solve the system f(y)=0. If an iterative linear
%   solver was specified (see KINSetOptions), KINSol uses a nonlinear Krylov
%   subspace projection method. The Newton-Krylov iterations are stopped
%   if either of the following conditions is satisfied:
%
%       ||f(y)||_L-infinity <= 0.01*fnormtol
%
%       ||y[i+1] - y[i]||_L-infinity <= scsteptol
%
%   However, if the current iterate satisfies the second stopping
%   criterion, it doesn't necessarily mean an approximate solution
%   has been found since the algorithm may have stalled, or the
%   user-specified step tolerance may be too large.
%
%   STRATEGY specifies the global strategy applied to the Newton step if it is
%   unsatisfactory. Valid choices are 'None' or 'LineSearch'.
%   YSCALE is a vector containing diagonal elements of scaling matrix for vector 
%   Y chosen so that the components of YSCALE*Y (as a matrix multiplication) all 
%   have about the same magnitude when Y is close to a root of f(y)
%   FSCALE is a vector containing diagonal elements of scaling matrix for f(y) 
%   chosen so that the components of FSCALE*f(Y) (as a matrix multiplication) 
%   all have roughly the same magnitude when u is not too near a root of f(y)
%
%   On return, status is one of the following:
%     0: KINSol succeeded
%     1: The initial y0 already satisfies the stopping criterion given above
%     2: Stopping tolerance on scaled step length satisfied
%    -1: An error occurred (see printed error message)
%
%   See also KINSetOptions, KINGetstats

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
% $Revision$Date: 2006/01/06 19:00:02 $

mode = 2;

[status, y] = kim(mode, y0, strategy, yscale, fscale);

