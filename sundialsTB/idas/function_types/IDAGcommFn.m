%IDAGcommFn - type for user provided communication function (BBDPre).
%
%IVP Problem
%
%   The function GCOMFUN must be defined as 
%        FUNCTION FLAG = GCOMFUN(T, YY, YP)
%   and can be used to perform all interprocess communication necessary
%   to evaluate the approximate residual function for the BBDPre
%   preconditioner module.
%   If a user data structure DATA was specified in IDAMalloc, then
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
%Adjoint Problem
%
%   The function GCOMFUNB must be defined either as
%        FUNCTION FLAG = GCOMFUNB(T, YY, YP, YYB, YPB)
%   or as
%        FUNCTION [FLAG, NEW_DATA] = GCOMFUNB(T, YY, YP, YYB, YPB, DATA)
%   depending on whether a user data structure DATA was specified in
%   IDAMalloc. 
%
%   The function GCOMFUNB must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also IDAGlocalFn, IDASetOptions
%
%   NOTES:
%     GCOMFUN and GCOMFUNB are specified through the GcommFn property in
%     IDASetOptions and are used only if the property PrecModule is set 
%     to 'BBDPre'.
%
%     Each call to GCOMFUN is preceded by a call to the residual function
%     DAEFUN with the same arguments T, YY, and YP (and YYB and YPB in the 
%     case of GCOMFUNB). Thus GCOMFUN can omit any communication done by 
%     DAEFUN if relevant to the evaluation of G by GLOCFUN. 
%     If all necessary communication was done by DAEFUN, GCOMFUN need 
%     not be provided.     

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/07/17 16:49:49 $
