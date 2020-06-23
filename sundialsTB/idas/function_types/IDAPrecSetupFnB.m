%IDAPrecSetupFnB - type for preconditioner setup function for backward problems.
%
%   The function PSETFUNB must be defined either as
%        FUNCTION FLAG = PSETFUNB(T,YY,YP,YYB,YPB,RRB,CJB)
%   or as
%        FUNCTION [FLAG,NEW_DATA] = PSETFUNB(T,YY,YP,YYB,YPB,RRB,CJB,DATA)
%   depending on whether a user data structure DATA was specified in
%   IDASetUserData.
%
%   See also IDAPrecSolveFnB, IDAPrecSetupFn, IDASetOptions
%
%   NOTE: PSETFUN and PSETFUNB are specified through the property
%   PrecSetupFn to IDASetOptions and are used only if the property
%   LinearSolver was set to 'GMRES', 'BiCGStab', or 'TFQMR'.

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
% $Revision$Date: 2007/08/21 17:38:44 $
