function status = IDAQuadReInit(yQ0, options)
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

mode = 12;

if nargin < 1
  error('Too few input arguments');
end

if nargin < 2
  options = [];
end
  
status = idm(mode, yQ0, options);
