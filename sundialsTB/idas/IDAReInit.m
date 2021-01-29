function status = IDAReInit(t0,yy0,yp0,options)
%IDAReInit reinitializes memory for IDAS.
%   where a prior call to IDAInit has been made with the same
%   problem size N. IDAReInit performs the same input checking
%   and initializations that IDAInit does, but it does no 
%   memory allocation, assuming that the existing internal memory 
%   is sufficient for the new problem.
%
%   Usage: IDAReInit ( T0, YY0, YP0 [, OPTIONS ] ) 
%
%   T0       is the initial value of t.
%   YY0      is the initial condition vector y(t0).  
%   YP0      is the initial condition vector y'(t0).  
%   OPTIONS  is an (optional) set of integration options, created with
%            the IDASetOptions function. 
%
%  See also: IDASetOptions, IDAInit

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
% $Revision$Date: 2007/08/21 17:38:42 $

mode = 11;

if nargin < 3
  error('Too few input arguments');
end

if nargin < 4
  options = [];
end

status = idm(mode, t0, yy0, yp0, options);
