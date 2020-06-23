function status = IDASet(varargin)
%IDASet changes optional input values during the integration.
%
%   Usage: IDASet('NAME1',VALUE1,'NAME2',VALUE2,...)
%
%   IDASet can be used to change some of the optional inputs during
%   the integration, i.e., without need for a solver reinitialization.
%   The property names accepted by IDASet are a subset of those valid
%   for IDASetOptions. Any unspecified properties are left unchanged.
%
%   IDASet with no input arguments displays all property names. 
%
%IDASet properties
%(See also the IDAS User Guide)
%
%UserData - problem data passed unmodified to all user functions.
%  Set VALUE to be the new user data.  
%RelTol - Relative tolerance
%  Set VALUE to the new relative tolerance
%AbsTol - absolute tolerance
%  Set VALUE to be either the new scalar absolute tolerance or
%  a vector of absolute tolerances, one for each solution component.
%StopTime - Stopping time
%  Set VALUE to be a new value for the independent variable past which
%  the solution is not to proceed. 

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
% $Revision$Date: 2007/08/21 17:38:43 $

if (nargin == 0)
  fprintf('        UserData\n');
  fprintf('\n');
  fprintf('          RelTol\n');
  fprintf('          AbsTol\n');
  fprintf('        StopTime\n');  
  fprintf('\n');
  return;
end

KeyNames = {
    'UserData'
    'RelTol'
    'AbsTol'
    'StopTime'
           };

options = idm_options(KeyNames,varargin{:});

mode = 33;

status = idm(mode, options);
