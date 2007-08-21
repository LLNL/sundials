function CVodeQuadInitB(idxB, fctQB, yQB0, optionsB)
%CVodeQuadInitB allocates and initializes memory for backward quadrature integration.
%
%   Usage: CVodeQuadInitB ( IDXB, QBFUN, YQB0 [, OPTIONS ] ) 
%
%   IDXB     is the index of the backward problem, returned by
%            CVodeInitB.
%   QBFUN    is a function defining the righ-hand sides of the
%            backward ODEs yQB' = fQB(t,y,yB).
%   YQB0     is the final conditions vector yQB(tB0).
%   OPTIONS  is an (optional) set of QUAD options, created with
%            the CVodeSetQuadOptions function. 
%
%   See also: CVodeInitB, CVodeSetQuadOptions, CVQuadRhsFnB 
%

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.2 $Date: 2007/05/11 18:51:32 $

mode = 6;

if nargin < 3
  error('Too few input arguments');
end

if nargin < 4
  optionsB = [];
end

idxB = idxB-1;
cvm(mode, idxB, fctQB, yQB0, optionsB);
