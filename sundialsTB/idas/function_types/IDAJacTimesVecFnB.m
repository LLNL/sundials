%IDAJacTimesVecFn - type for Jacobian times vector function for backward problems.
%
%   The function JTVFUNB must be defined either as
%        FUNCTION [JVB,FLAG] = JTVFUNB(T,YY,YP,YYB,YPB,RRB,VB,CJB)
%   or as
%        FUNCTION [JVB,FLAG,NEW_DATA] = JTVFUNB(T,YY,YP,YYB,YPB,RRB,VB,CJB,DATA)
%   depending on whether a user data structure DATA was specified in
%   IDAInit. In either case, it must return the vector JVB, the
%   product of the Jacobian (dfB/dyyB + cj * dfB/dypB) and a vector
%   vB. The input argument RRB contains the current value of f(t,yy,yp,yyB,ypB).
%
%   The function JTVFUNB must set FLAG=0 if successful, or FLAG~=0 if
%   a failure occurred.
%
%   See also IDASetOptions
%
%   NOTE: JTVFUNB is specified through the property JacobianFn to IDASetOptions
%   and is used only if the property LinearSolver was set to 'GMRES', 'BiCGStab',
%   or 'TFQMR'.

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
