function CVodeReInitB(idxB, tB0, yB0, optionsB)
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
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.3 $Date: 2007/05/11 18:51:32 $

mode = 15;

if nargin < 3
  error('Too few input arguments');
end

if nargin < 4
  optionsB = [];
end

idxB = idxB-1;
cvm(mode,idxB,tB0,yB0,optionsB);
