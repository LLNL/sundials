function status = IDAQuadInitB(idxB, fctQB, yQB0, optionsB)
%IDAQuadInitB allocates and initializes memory for backward quadrature integration.
%
%   Usage: IDAQuadInitB ( IDXB, QBFUN, YQB0 [, OPTIONS ] ) 
%
%   IDXB     is the index of the backward problem, returned by
%            IDAInitB.
%   QBFUN    is a function defining the righ-hand sides of the
%            backward ODEs yQB' = fQB(t,y,yB).
%   YQB0     is the final conditions vector yQB(tB0).
%   OPTIONS  is an (optional) set of QUAD options, created with
%            the IDASetQuadOptions function. 
%
%   See also: IDAInitB, IDASetQuadOptions, IDAQuadRhsFnB 
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

mode = 6;

if nargin < 3
  error('Too few input arguments');
end

if nargin < 4
  optionsB = [];
end

idxB = idxB-1;
status = idm(mode, idxB, fctQB, yQB0, optionsB);
