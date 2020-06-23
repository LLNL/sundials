%CVDenseJacFnB - type for user provided dense Jacobian function for backward problems.
%
%   The function DJACFUNB must be defined either as
%        FUNCTION [JB, FLAG] = DJACFUNB(T, Y, YB, FYB)
%   or as
%        FUNCTION [JB, FLAG, NEW_DATA] = DJACFUNB(T, Y, YB, FYB, DATA)
%   depending on whether a user data structure DATA was specified in
%   CVodeInit. In either case, it must return the matrix JB, the
%   Jacobian of fB(t,y,yB), with respect to yB. The input argument
%   FYB contains the current value of f(t,y,yB).
%
%   The function DJACFUNB must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also CVodeSetOptions
%
%   NOTE: DJACFUNB is specified through the property JacobianFn to
%   CVodeSetOptions and is used only if the property LinearSolver
%   was set to 'Dense'.

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
% $Revision$Date: 2007/05/11 18:51:33 $
