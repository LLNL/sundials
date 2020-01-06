%CVGlocalFnB - type for user provided RHS approximation function (BBDPre) for backward problems.
%
%   The function GLOCFUNB must be defined either as
%        FUNCTION [GLOCB, FLAG] = GLOCFUNB(T,Y,YB)
%   or as
%        FUNCTION [GLOCB, FLAG, NEW_DATA] = GLOCFUNB(T,Y,YB,DATA)
%   depending on whether a user data structure DATA was specified in
%   CVodeInit. In either case, it must return the vector GLOCB
%   corresponding to an approximation to fB(t,y,yB).
%
%   The function GLOCFUNB must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also CVGcommFnB, CVodeSetOptions
%
%   NOTE: GLOCFUNB is specified through the GlocalFn property in CVodeSetOptions
%   and is used only if the property PrecModule is set to 'BBDPre'.

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
