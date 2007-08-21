function IDAQuadReInit(yQ0, options)
%IDAQuadReInit reinitializes IDAS's quadrature-related memory
%   assuming it has already been allocated in prior calls to IDAInit 
%   and IDAQuadInit.
%
%   Usage: IDAQuadReInit ( YQ0 [, OPTIONS ] ) 
%
%   YQ0      Initial conditions for quadrature variables yQ(t0).
%   OPTIONS  is an (optional) set of QUAD options, created with
%            the IDASetQuadOptions function. 
%
%   See also: IDASetQuadOptions, IDAQuadInit

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.1 $Date: 2007/05/11 18:51:32 $

mode = 12;

if nargin < 1
  error('Too few input arguments');
end

if nargin < 2
  options = [];
end
  
idm(mode, yQ0, options);
