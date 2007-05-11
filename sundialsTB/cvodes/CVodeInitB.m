function idxB = CVodeInitB(fctB, tB0, yB0, optionsB)
%CVodeInitB allocates and initializes backward memory for CVODES.
%
%   Usage:   IDXB = CVodeInitB ( FCTB, TB0, YB0 [, OPTIONSB] )
%
%   FCTB     is a function defining the adjoint ODE right-hand side.
%            This function must return a vector containing the current 
%            value of the adjoint ODE righ-hand side.
%   TB0      is the final value of t.
%   YB0      is the final condition vector yB(tB0).  
%   OPTIONSB is an (optional) set of integration options, created with
%            the CVodeSetOptions function. 
%
%   CVodeInitB returns the index IDXB associated with this backward
%   problem. This index must be passed as an argument to any subsequent
%   functions related to this backward problem.
%
%   See also: CVRhsFnB
%

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/07/07 19:08:40 $

mode = 5;

if nargin < 3
  error('Too few input arguments');
end

if nargin < 4
  optionsB = [];
end

idxB = cvm(mode, fctB, tB0, yB0, optionsB);
idxB = idxB+1;