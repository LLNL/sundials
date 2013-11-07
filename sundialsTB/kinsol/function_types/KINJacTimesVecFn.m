%KINJacTimesVecFn - type for user provided Jacobian times vector function.
%
%   The function JTVFUN must be defined as 
%        FUNCTION [JV, NEW_Y, FLAG] = JTVFUN(Y, V, NEW_Y)
%   and must return a vector JV corresponding to the product of the 
%   Jacobian of f(y) with the vector v. On input, NEW_Y indicates if
%   the iterate has been updated in the interim. JV must be update
%   or reevaluated, if appropriate, unless NEW_Y=false. This flag must
%   be reset by the user.
%   If a user data structure DATA was specified in KINInit, then
%   JTVFUN must be defined as
%        FUNCTION [JV, NEW_Y, FLAG, NEW_DATA] = JTVFUN(Y, V, NEW_Y, DATA)
%   If the local modifications to the user data structure are needed in
%   other user-provided functions then, besides setting the vector JV, and
%   flags NEW_Y and FLAG, the JTVFUN function must also set NEW_DATA. Otherwise, 
%   it should set NEW_DATA=[] (do not set NEW_DATA = DATA as it would lead to
%   unnecessary copying).
%
%   If successful, FLAG should be set to 0. If an error occurs, FLAG should
%   be set to a nonzero value.
%
%   See also KINSetOptions
%
%   NOTE: JTVFUN is specified through the property JacobianFn to KINSetOptions 
%   and is used only if the property LinearSolver was set to 'GMRES' or 'BiCGStab'.

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
% $Revision: 1.2 $Date: 2007/05/11 18:48:46 $
