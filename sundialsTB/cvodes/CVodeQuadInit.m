function CVodeQuadInit(fctQ, yQ0, options)
%CVodeQuadInit allocates and initializes memory for quadrature integration.
%
%   Usage: CVodeQuadInit ( QFUN, YS0 [, OPTIONS ] ) 
%
%   QFUN     is a function defining the righ-hand sides of the quadrature
%            ODEs yQ' = fQ(t,y).
%   YQ0      is the initial conditions vector yQ(t0).
%   OPTIONS  is an (optional) set of QUAD options, created with
%            the CVodeSetQuadOptions function. 
%
%   See also: CVodeSetQuadOptions, CVQuadRhsFn 

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/03/07 01:19:50 $

mode = 2;

if nargin < 2
  error('Too few input arguments');
end

if nargin < 3
  options = [];
end

cvm(mode, fctQ, yQ0, options);
