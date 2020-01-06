%IDAJacTimesVecFn - type for Jacobian times vector function.
%
%   The function JTVFUN must be defined as 
%        FUNCTION [JV, FLAG] = JTVFUN(T,YY,YP,RR,V,CJ)
%   and must return a vector JV corresponding to the product of the 
%   Jacobian ( df/dyy + cj * df/dyp ) with the vector v.
%   The input argument RR contains the current value of f(t,yy,yp).
%   If a user data structure DATA was specified in IDAInit, then
%   JTVFUN must be defined as
%        FUNCTION [JV, FLAG, NEW_DATA] = JTVFUN(T,YY,YP,RR,V,CJ,DATA)
%   If the local modifications to the user data structure are needed in
%   other user-provided functions then, besides setting the vector JV,
%   the JTVFUN function must also set NEW_DATA. Otherwise, it should set
%   NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead to
%   unnecessary copying).
%
%   The function JTVFUN must set FLAG=0 if successful, or FLAG~=0 if
%   a failure occurred.
%
%   See also IDASetOptions
%
%   NOTE: JTVFUN is specified through the property JacobianFn to 
%   IDASetOptions and is used only if the property LinearSolver 
%   was set to 'GMRES', 'BiCGStab', or 'TFQMR'.

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
