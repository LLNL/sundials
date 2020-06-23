%CVQuadRhsFn - type for user provided quadrature RHS function.
%
%   The function ODEQFUN must be defined as 
%        FUNCTION [YQD, FLAG] = ODEQFUN(T,Y)
%   and must return a vector YQD corresponding to fQ(t,y), the integrand
%   for the integral to be evaluated.
%   If a user data structure DATA was specified in CVodeInit, then
%   ODEQFUN must be defined as
%        FUNCTION [YQD, FLAG, NEW_DATA] = ODEQFUN(T,Y,DATA)
%   If the local modifications to the user data structure are needed in
%   other user-provided functions then, besides setting the vector YQD,
%   the ODEQFUN function must also set NEW_DATA. Otherwise, it should set
%   NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead to
%   unnecessary copying).
%
%   The function ODEQFUN must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also CVodeQuadInit

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
