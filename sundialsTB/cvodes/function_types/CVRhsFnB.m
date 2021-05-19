%CVRhsFnB - type for user provided RHS function for backward problems.
%
%   The function ODEFUNB must be defined either as
%        FUNCTION [YBD, FLAG] = ODEFUNB(T,Y,YB)
%   or as
%        FUNCTION [YBD, FLAG, NEW_DATA] = ODEFUNB(T,Y,YB,DATA)
%   depending on whether a user data structure DATA was specified in
%   CVodeInit. In either case, it must return the vector YBD
%   corresponding to fB(t,y,yB).
%
%   The function ODEFUNB must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also CVodeInitB
%

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
