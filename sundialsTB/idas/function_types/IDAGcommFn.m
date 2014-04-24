%IDAGcommFn - type for communication function (BBDPre).
%
%   The function GCOMFUN must be defined as 
%        FUNCTION FLAG = GCOMFUN(T, YY, YP)
%   and can be used to perform all interprocess communication necessary
%   to evaluate the approximate residual function for the BBDPre
%   preconditioner module.
%   If a user data structure DATA was specified in IDAInit, then
%   GCOMFUN must be defined as
%        FUNCTION [FLAG, NEW_DATA] = GCOMFUN(T, YY, YP, DATA)
%   If the local modifications to the user data structure are needed 
%   in other user-provided functions then the GCOMFUN function must also 
%   set NEW_DATA. Otherwise, it should set NEW_DATA=[] (do not set 
%   NEW_DATA = DATA as it would lead to unnecessary copying).
%
%   The function GCOMFUN must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also IDAGlocalFn, IDASetOptions
%
%   NOTES:
%     GCOMFUN is specified through the GcommFn property in IDASetOptions 
%     and is used only if the property PrecModule is set to 'BBDPre'.
%
%     Each call to GCOMFUN is preceded by a call to the residual function
%     DAEFUN with the same arguments T, YY, and YP.
%     Thus GCOMFUN can omit any communication done by DAEFUN if relevant 
%     to the evaluation of G by GLOCFUN. If all necessary communication 
%     was done by DAEFUN, GCOMFUN need not be provided.     

% Radu Serban <radu@llnl.gov>
% LLNS Copyright Start
% Copyright (c) 2014, Lawrence Livermore National Security
% This work was performed under the auspices of the U.S. Department 
% of Energy by Lawrence Livermore National Laboratory in part under 
% Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
% Produced at the Lawrence Livermore National Laboratory.
% All rights reserved.
% For details, see the LICENSE file.
% LLNS Copyright End
% $Revision$Date: 2007/08/21 17:38:44 $
