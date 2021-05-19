%CVPrecSolveFn - type for user provided preconditioner solve function.
%
%   The user-supplied preconditioner solve function PSOLFN
%   is to solve a linear system P z = r in which the matrix P is
%   one of the preconditioner matrices P1 or P2, depending on the
%   type of preconditioning chosen.
%
%   The function PSOLFUN must be defined as 
%        FUNCTION [Z, FLAG] = PSOLFUN(T,Y,FY,R)
%   and must return a vector Z containing the solution of Pz=r.
%   If PSOLFUN was successful, it must return FLAG=0. For a recoverable 
%   error (in which case the step will be retried) it must set FLAG to a 
%   positive value. If an unrecoverable error occurs, it must set FLAG
%   to a negative value, in which case the integration will be halted.
%   The input argument FY contains the current value of f(t,y).
%
%   If a user data structure DATA was specified in CVodeInit, then
%   PSOLFUN must be defined as
%        FUNCTION [Z, FLAG, NEW_DATA] = PSOLFUN(T,Y,FY,R,DATA)
%   If the local modifications to the user data structure are needed in
%   other user-provided functions then, besides setting the vector Z and
%   the flag FLAG, the PSOLFUN function must also set NEW_DATA. Otherwise,
%   it should set NEW_DATA=[] (do not set NEW_DATA = DATA as it would
%   lead to unnecessary copying).
%
%   See also CVPrecSetupFn, CVodeSetOptions
%
%   NOTE: PSOLFUN is specified through the property PrecSolveFn to 
%   CVodeSetOptions and is used only if the property LinearSolver was
%   set to 'GMRES', 'BiCGStab', or 'TFQMR' and if the property PrecType
%   is not 'None'.

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
% $Revision$Date: 2007/05/11 18:51:33 $
