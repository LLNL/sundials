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
% SUNDIALS Copyright Start
% Copyright (c) 2002-2020, Lawrence Livermore National Security
% and Southern Methodist University.
% All rights reserved.
%
% See the top-level LICENSE and NOTICE files for details.
%
% SPDX-License-Identifier: BSD-3-Clause
% SUNDIALS Copyright End
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
