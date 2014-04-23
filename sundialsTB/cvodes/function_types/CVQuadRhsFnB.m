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
% LLNS Start Copyright
% Copyright (c) 2013, Lawrence Livermore National Security
% This work was performed under the auspices of the U.S. Department 
% of Energy by Lawrence Livermore National Laboratory in part under 
% Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
% Produced at the Lawrence Livermore National Laboratory.
% All rights reserved.
% For details, see the LICENSE file.
% LLNS End Copyright
% $Revision$Date: 2007/05/11 18:51:33 $
