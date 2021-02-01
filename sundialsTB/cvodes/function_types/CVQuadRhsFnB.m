%CVQuadRhsFnB - type for user provided quadrature RHS function for backward problems
%
%   The function ODEQFUNB must be defined either as
%        FUNCTION [YQBD, FLAG] = ODEQFUNB(T,Y,YB)
%   or as
%        FUNCTION [YQBD, FLAG, NEW_DATA] = ODEQFUNB(T,Y,YB,DATA)
%   depending on whether a user data structure DATA was specified in
%   CVodeInit. In either case, it must return the vector YQBD
%   corresponding to fQB(t,y,yB), the integrand for the integral to be 
%   evaluated on the backward phase.
%
%   The function ODEQFUNB must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also CVodeQuadInitB

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
