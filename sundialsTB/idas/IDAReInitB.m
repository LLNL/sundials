function status = IDAReInitB(idxB,tB0,yyB0,ypB0,optionsB)
%IDAReInitB allocates and initializes backward memory for IDAS.
%   where a prior call to IDAInitB has been made with the same
%   problem size NB. IDAReInitB performs the same input checking
%   and initializations that IDAInitB does, but it does no 
%   memory allocation, assuming that the existing internal memory 
%   is sufficient for the new problem.
%
%   Usage:   IDAReInitB ( IDXB, TB0, YYB0, YPB0 [, OPTIONSB] )
%
%   IDXB     is the index of the backward problem, returned by
%            IDAInitB.
%   TB0      is the final value of t.
%   YYB0     is the final condition vector yB(tB0).  
%   YPB0     is the final condition vector yB'(tB0).
%   OPTIONSB is an (optional) set of integration options, created with
%            the IDASetOptions function. 
%
%   See also: IDASetOptions, IDAInitB
%

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
% $Revision$Date: 2007/08/21 17:38:42 $

mode = 15;

if nargin < 4
  error('Too few input arguments');
end

if nargin < 5
  optionsB = [];
end

idxB = idxB-1;
status = idm(mode, idxB, tB0, yyB0, ypB0, optionsB);
