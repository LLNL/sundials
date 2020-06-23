function status = CVodeQuadReInitB(idxB, yQB0, optionsB)
%CVodeQuadReInitB reinitializes memory for backward quadrature integration.
%
%   Usage: CVodeQuadReInitB ( IDXB, YS0 [, OPTIONS ] ) 
%
%   IDXB     is the index of the backward problem, returned by
%            CVodeInitB.
%   YQB0     is the final conditions vector yQB(tB0).
%   OPTIONS  is an (optional) set of QUAD options, created with
%            the CVodeSetQuadOptions function. 
%
%   See also: CVodeSetQuadOptions, CVodeReInitB, CVodeQuadInitB
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
% $Revision$Date: 2007/05/11 18:51:32 $

mode = 16;

if nargin < 2
  error('Too few input arguments');
end

if nargin < 3
  optionsB = [];
end
  
idxB = idxB-1;
status = cvm(mode, idxB, yQB0, optionsB);
