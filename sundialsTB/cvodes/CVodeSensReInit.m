function status = CVodeSensReInit(yS0, options)
%CVodeSensReInit reinitializes CVODES's FSA-related memory
%   assuming it has already been allocated in prior calls to CVodeInit 
%   and CVodeSensInit.
%   The number of sensitivities Ns is assumed to be unchanged since the 
%   previous call to CVodeSensInit.
%
%   Usage: CVodeSensReInit ( YS0 [, OPTIONS ] ) 
%
%   YS0      Initial conditions for sensitivity variables.
%            YS0 must be a matrix with N rows and Ns columns, where N is the problem
%            dimension and Ns the number of sensitivity systems.
%   OPTIONS  is an (optional) set of FSA options, created with
%            the CVodeSensSetOptions function. 
%
%   See also: CVodeSensSetOptions, CVodeReInit, CVodeSensInit

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

mode = 13;

if nargin < 1
  error('Too few input arguments');
end

if nargin < 2
  options = [];
end

status = cvm(mode, yS0, options);
