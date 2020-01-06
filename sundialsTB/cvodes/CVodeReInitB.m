function status = CVodeReInitB(idxB, tB0, yB0, optionsB)
%CVodeReInitB re-initializes backward memory for CVODES.
%   where a prior call to CVodeInitB has been made with the same
%   problem size NB. CVodeReInitB performs the same input checking
%   and initializations that CVodeInitB does, but it does no 
%   memory allocation, assuming that the existing internal memory 
%   is sufficient for the new problem.
%
%   Usage:   CVodeReInitB ( IDXB, TB0, YB0 [, OPTIONSB] )
%
%   IDXB     is the index of the backward problem, returned by
%            CVodeInitB.
%   TB0      is the final value of t.
%   YB0      is the final condition vector yB(tB0).  
%   OPTIONSB is an (optional) set of integration options, created with
%            the CVodeSetOptions function. 
%
%   See also: CVodeSetOptions, CVodeInitB
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
% $Revision$Date: 2007/08/21 17:42:38 $

mode = 15;

if nargin < 3
  error('Too few input arguments');
end

if nargin < 4
  optionsB = [];
end

idxB = idxB-1;
status = cvm(mode,idxB,tB0,yB0,optionsB);
