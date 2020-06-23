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

mode = 2;

if nargin < 2
  error('Too few input arguments');
end

if nargin < 3
  options = [];
end

status = idm(mode, fctQ, yQ0, options);
