function status = CVodeQuadReInit(yQ0, options)
%CVodeQuadReInit reinitializes CVODES's quadrature-related memory
%   assuming it has already been allocated in prior calls to CVodeInit 
%   and CVodeQuadInit.
%
%   Usage: CVodeQuadReInit ( YQ0 [, OPTIONS ] ) 
%
%   YQ0      Initial conditions for quadrature variables yQ(t0).
%   OPTIONS  is an (optional) set of QUAD options, created with
%            the CVodeSetQuadOptions function. 
%
%   See also: CVodeSetQuadOptions, CVodeQuadInit

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
% $Revision$Date: 2007/05/11 18:51:32 $

mode = 12;

if nargin < 1
  error('Too few input arguments');
end

if nargin < 2
  options = [];
end
  
status = cvm(mode, yQ0, options);
