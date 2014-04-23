function status = IDAQuadInit(fctQ, yQ0, options)
%IDAQuadInit allocates and initializes memory for quadrature integration.
%
%   Usage: IDAQuadInit ( QFUN, YQ0 [, OPTIONS ] ) 
%
%   QFUN     is a function defining the righ-hand sides of the quadrature
%            ODEs yQ' = fQ(t,y).
%   YQ0      is the initial conditions vector yQ(t0).
%   OPTIONS  is an (optional) set of QUAD options, created with
%            the IDASetQuadOptions function. 
%
%   See also: IDASetQuadOptions, IDAQuadRhsFn 

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

mode = 2;

if nargin < 2
  error('Too few input arguments');
end

if nargin < 3
  options = [];
end

status = idm(mode, fctQ, yQ0, options);
