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
% LLNS Copyright Start
% Copyright (c) 2014, Lawrence Livermore National Security
% This work was performed under the auspices of the U.S. Department 
% of Energy by Lawrence Livermore National Laboratory in part under 
% Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
% Produced at the Lawrence Livermore National Laboratory.
% All rights reserved.
% For details, see the LICENSE file.
% LLNS Copyright End
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
