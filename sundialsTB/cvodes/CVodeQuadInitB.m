function status = CVodeQuadInitB(idxB, fctQB, yQB0, optionsB)
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
% SUNDIALS Copyright Start
% Copyright (c) 2002-2021, Lawrence Livermore National Security
% and Southern Methodist University.
% All rights reserved.
%
% See the top-level LICENSE and NOTICE files for details.
%
% SPDX-License-Identifier: BSD-3-Clause
% SUNDIALS Copyright End
% $Revision$Date: 2007/08/21 17:42:38 $

mode = 6;

if nargin < 3
  error('Too few input arguments');
end

if nargin < 4
  optionsB = [];
end

idxB = idxB-1;
status = cvm(mode, idxB, fctQB, yQB0, optionsB);
