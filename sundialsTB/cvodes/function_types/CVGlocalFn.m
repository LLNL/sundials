%CVGlocalFn - type for user provided RHS approximation function (BBDPre).
%
%   The function GLOCFUN must be defined as 
%        FUNCTION [GLOC, FLAG] = GLOCFUN(T,Y)
%   and must return a vector GLOC corresponding to an approximation to f(t,y)
%   which will be used in the BBDPRE preconditioner module. The case where
%   G is mathematically identical to F is allowed.
%   If a user data structure DATA was specified in CVodeInit, then
%   GLOCFUN must be defined as
%        FUNCTION [GLOC, FLAG, NEW_DATA] = GLOCFUN(T,Y,DATA)
%   If the local modifications to the user data structure are needed 
%   in other user-provided functions then, besides setting the vector G,
%   the GLOCFUN function must also set NEW_DATA. Otherwise, it should set
%   NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead to
%   unnecessary copying).
%
%   The function GLOCFUN must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also CVGcommFn, CVodeSetOptions
%
%   NOTE: GLOCFUN is specified through the GlocalFn property in CVodeSetOptions
%   and is used only if the property PrecModule is set to 'BBDPre'.

% Radu Serban <radu@llnl.gov>
% LLNS Start Copyright
% Copyright (c) 2013, Lawrence Livermore National Security
% This work was performed under the auspices of the U.S. Department 
% of Energy by Lawrence Livermore National Laboratory in part under 
% Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
% Produced at the Lawrence Livermore National Laboratory.
% All rights reserved.
% For details, see the LICENSE file.
% LLNS End Copyright
% $Revision: 1.2 $Date: 2007/05/11 18:51:33 $
