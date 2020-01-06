function options = CVodeQuadSetOptions(varargin)
%CVodeQuadSetOptions creates an options structure for quadrature integration with CVODES.
%
%   Usage: OPTIONS = CVodeQuadSetOptions('NAME1',VALUE1,'NAME2',VALUE2,...)
%          OPTIONS = CVodeQuadSetOptions(OLDOPTIONS,'NAME1',VALUE1,...)
%
%   OPTIONS = CVodeQuadSetOptions('NAME1',VALUE1,'NAME2',VALUE2,...) creates 
%   a CVODES options structure OPTIONS in which the named properties have 
%   the specified values. Any unspecified properties have default values. 
%   It is sufficient to type only the leading characters that uniquely 
%   identify the property. Case is ignored for property names. 
%   
%   OPTIONS = CVodeQuadSetOptions(OLDOPTIONS,'NAME1',VALUE1,...) alters an 
%   existing options structure OLDOPTIONS.
%   
%   CVodeQuadSetOptions with no input arguments displays all property names 
%   and their possible values.
%   
%CVodeQuadSetOptions properties
%(See also the CVODES User Guide)
%
%ErrControl - Error control strategy for quadrature variables [ {false} | true ]
%   Specifies whether quadrature variables are included in the error test.
%RelTol - Relative tolerance for quadrature variables [ scalar {1e-4} ]
%   Specifies the relative tolerance for quadrature variables. This parameter is
%   used only if ErrControl = true.
%AbsTol - Absolute tolerance for quadrature variables [ scalar or vector {1e-6} ]
%   Specifies the absolute tolerance for quadrature variables. This parameter is
%   used only if ErrControl = true.
%
%SensDependent - Backward problem depending on sensitivities [ {false} | true ]
%   Specifies whether the backward problem quadrature right-hand side depends
%   on forward sensitivites. If SUNTRUE, the right-hand side function provided for
%   this backward problem must have the appropriate type (see CVQuadRhsFnB).
%
%
%   See also
%        CVodeQuadInit, CVodeQuadReInit.
%        CVodeQuadInitB, CVodeQuadReInitB

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
% $Revision$Date: 2007/05/11 21:42:52 $

% If called without input and output arguments, print out the possible keywords

if (nargin == 0) && (nargout == 0)
  fprintf('      ErrControl: [ {false} | true ]\n');
  fprintf('          RelTol: [ positive scalar {1e-4} ]\n');
  fprintf('          AbsTol: [ positive scalar or vector {1e-6} ]\n');
  fprintf('\n');
  fprintf('   SensDependent: [ {false} | true ]\n');
  fprintf('\n');
  return;
end

KeyNames = {
    'ErrControl'
    'RelTol'
    'AbsTol'
    'SensDependent'
    };

options = cvm_options(KeyNames,varargin{:});
