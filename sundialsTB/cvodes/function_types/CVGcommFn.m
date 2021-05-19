%CVGcommFn - type for user provided communication function (BBDPre).
%
%   The function GCOMFUN must be defined as 
%        FUNCTION FLAG = GCOMFUN(T, Y)
%   and can be used to perform all interprocess communication necessary
%   to evaluate the approximate right-hand side function for the BBDPre
%   preconditioner module.
%   If a user data structure DATA was specified in CVodeInit, then
%   GCOMFUN must be defined as
%        FUNCTION [FLAG, NEW_DATA] = GCOMFUN(T, Y, DATA)
%   If the local modifications to the user data structure are needed 
%   in other user-provided functions then the GCOMFUN function must also 
%   set NEW_DATA. Otherwise, it should set NEW_DATA=[] (do not set 
%   NEW_DATA = DATA as it would lead to unnecessary copying).
%
%   The function GCOMFUN must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also CVGlocalFn, CVodeSetOptions
%
%   NOTES:
%     GCOMFUN is specified through the GcommFn property in CVodeSetOptions
%     and is used only if the property PrecModule is set to 'BBDPre'.
%
%     Each call to GCOMFUN is preceded by a call to the RHS function
%     ODEFUN with the same arguments T and Y. Thus GCOMFUN can omit
%     any communication done by ODEFUN if relevant to the evaluation
%     of G by GLOCFUN. If all necessary communication was done by ODEFUN,
%     GCOMFUN need not be provided.     

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
