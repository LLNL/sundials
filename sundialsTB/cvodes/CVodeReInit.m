function status = CVodeReInit(t0, y0, options)
%CVodeReInit reinitializes memory for CVODES
%   where a prior call to CVodeInit has been made with the same
%   problem size N. CVodeReInit performs the same input checking
%   and initializations that CVodeInit does, but it does no 
%   memory allocation, assuming that the existing internal memory 
%   is sufficient for the new problem.
%
%   Usage: CVodeReInit ( T0, Y0 [, OPTIONS ] ) 
%
%   T0       is the initial value of t.
%   Y0       is the initial condition vector y(t0).  
%   OPTIONS  is an (optional) set of integration options, created with
%            the CVodeSetOptions function. 
%
%   See also: CVodeSetOptions, CVodeInit

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
% $Revision$Date: 2007/05/16 17:12:56 $

mode = 11;

if nargin < 2
  error('Too few input arguments');
end

if nargin < 3
  options = [];
end

status = cvm(mode, t0, y0, options);
