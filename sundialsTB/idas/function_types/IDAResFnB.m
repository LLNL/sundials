%IDAResFnb - type for residual function for backward problems
%
%   The function DAEFUNB must be defined either as
%        FUNCTION [RB, FLAG] = DAEFUNB(T, YY, YP, YYB, YPB)
%   or as
%        FUNCTION [RB, FLAG, NEW_DATA] = DAEFUNB(T, YY, YP, YYB, YPB, DATA)
%   depending on whether a user data structure DATA was specified in
%   IDAInit. In either case, it must return the vector RB
%   corresponding to fB(t,yy,yp,yyB,ypB).
%
%   The function DAEFUNB must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also IDAInitB, IDARhsFn

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
% $Revision: 1.2 $Date: 2007/08/21 17:38:44 $
