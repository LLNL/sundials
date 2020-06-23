function [si, status] = KINGetStats()
%KINGetStats returns statistics for the main KINSOL solver and the linear
%solver used.
%
%   Usage: STATS = KINGetStats
%
%Fields in the structure STATS
%
%o nfe    - total number evaluations of the nonlinear system function SYSFUN
%o nni    - total number of nonlinear iterations
%o nbcf   - total number of beta-condition failures
%o nbops  - total number of backtrack operations (step length adjustments) 
%           performed by the line search algorithm 
%o fnorm  - scaled norm of the nonlinear system function f(y) evaluated at the
%           current iterate: ||fscale*f(y)||_L2
%o step   - scaled norm (or length) of the step used during the previous 
%           iteration: ||uscale*p||_L2
%o LSInfo - structure with linear solver statistics
%
%The structure LSinfo has different fields, depending on the linear solver used.
%
%  Fields in LSinfo for the 'Dense' linear solver
%
%o name - 'Dense'
%o njeD - number of Jacobian evaluations
%o nfeD - number of right-hand side function evaluations for difference-quotient
%         Jacobian approximation
%
%  Fields in LSinfo for the 'GMRES' or 'BiCGStab' linear solver
%
%o name - 'GMRES' or 'BiCGStab'
%o nli  - number of linear solver iterations
%o npe  - number of preconditioner setups
%o nps  - number of preconditioner solve function calls
%o ncfl - number of linear system convergence test failures
%

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
% $Revision$Date: 2006/01/06 19:00:02 $


mode = 3;
[si, status] = kim(mode);
