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
% LLNS Copyright Start
% Copyright (c) 2014, Lawrence Livermore National Security
% This work was performed under the auspices of the U.S. Department 
% of Energy by Lawrence Livermore National Laboratory in part under 
% Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
% Produced at the Lawrence Livermore National Laboratory.
% All rights reserved.
% For details, see the LICENSE file.
% LLNS Copyright End
% $Revision$Date: 2007/05/11 18:51:33 $
