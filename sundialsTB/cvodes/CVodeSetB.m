function status = CVodeSetB(idxB, varargin)
%CVodeSetB changes optional input values during the integration.
%
%   Usage: CVodeSetB( IDXB, 'NAME1',VALUE1,'NAME2',VALUE2,... )
%
%   CVodeSetB can be used to change some of the optional inputs for
%   the backward problem identified by IDXB during the backward
%   integration, i.e., without need for a solver reinitialization.
%   The property names accepted by CVodeSet are a subset of those valid
%   for CVodeSetOptions. Any unspecified properties are left unchanged.
%   
%   CVodeSetB with no input arguments displays all property names.
%
%CVodeSetB properties
%(See also the CVODES User Guide)
%
%UserData - problem data passed unmodified to all user functions.
%  Set VALUE to be the new user data.  
%RelTol - Relative tolerance
%  Set VALUE to the new relative tolerance
%AbsTol - absolute tolerance
%  Set VALUE to be either the new scalar absolute tolerance or
%  a vector of absolute tolerances, one for each solution component.

% Radu Serban <radu@llnl.gov>
% LLNS Copyright Start
% Copyright (c) 2014, Lawrence Livermore National Security
% This work was performed under the auspices of the U.S. Department 
% of Energy by Lawrence Livermore National Laboratory in part under 
% Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
% Produced at the Lawrence Livermore National Laboratory.
% All rights reserved.
% For details, see the LICENSE file.
% LLNS Copyright End
% $Revision$Date: 2007/05/16 17:12:56 $

if (nargin == 0)
  fprintf('        UserData\n');
  fprintf('\n');
  fprintf('          RelTol\n');
  fprintf('          AbsTol\n');
  fprintf('\n');
  return;
end

KeyNames = {
    'UserData'
    'RelTol'
    'AbsTol'
           };

options = cvm_options(KeyNames,varargin{:});

mode = 34;

status = cvm(mode, idxB, options);
