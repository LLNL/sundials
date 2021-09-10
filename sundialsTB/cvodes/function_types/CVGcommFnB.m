%CVGcommFn - type for user provided communication function (BBDPre) for backward problems.
%
%   The function GCOMFUNB must be defined either as
%        FUNCTION FLAG = GCOMFUNB(T, Y, YB)
%   or as
%        FUNCTION [FLAG, NEW_DATA] = GCOMFUNB(T, Y, YB, DATA)
%   depending on whether a user data structure DATA was specified in
%   CVodeInit. 
%
%   The function GCOMFUNB must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also CVGlocalFnB, CVodeSetOptions
%
%   NOTES:
%     GCOMFUNB is specified through the GcommFn property in CVodeSetOptions
%     and is used only if the property PrecModule is set to 'BBDPre'.
%
%     Each call to GCOMFUNB is preceded by a call to the RHS function
%     ODEFUNB with the same arguments T, Y, and YB. Thus GCOMFUNB can
%     omit any communication done by ODEFUNB if relevant to the evaluation
%     of G by GLOCFUNB. If all necessary communication was done by ODEFUNB,
%     GCOMFUNB need not be provided.     

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
