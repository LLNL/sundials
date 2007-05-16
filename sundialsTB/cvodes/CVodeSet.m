function CVodeSet(varargin)
%CVodeSet changes optional input values during the integration.
%
%   Usage: CVodeSet('NAME1',VALUE1,'NAME2',VALUE2,...)
%
%   CVodeSet can be used to change some of the optional inputs during
%   the integration, i.e., without need for a solver reinitialization.
%   The property names accepted by CVodeSet are a subset of those valid
%   for CVodeSetOptions. Any unspecified properties are left unchanged.
%
%   CVodeSet with no input arguments displays all property names. 
%
%CVodeSet properties
%(See also the CVODES User Guide)
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
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.1 $Date: 2007/05/11 18:51:31 $

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

options = cvm_options(KeyNames,varargin{:});

mode = 33;

cvm(mode, options);