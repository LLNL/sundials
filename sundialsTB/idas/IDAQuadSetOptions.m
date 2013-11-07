function options = IDAQuadSetOptions(varargin)
%IDAQuadSetOptions creates an options structure for IDAS.
%
%   Usage: OPTIONS = IDAQuadSetOptions('NAME1',VALUE1,'NAME2',VALUE2,...)
%          OPTIONS = IDAQuadSetOptions(OLDOPTIONS,'NAME1',VALUE1,...)
%
%   OPTIONS = IDAQuadSetOptions('NAME1',VALUE1,'NAME2',VALUE2,...) creates 
%   an IDAS options structure OPTIONS in which the named properties have 
%   the specified values. Any unspecified properties have default values. 
%   It is sufficient to type only the leading characters that uniquely 
%   identify the property. Case is ignored for property names. 
%   
%   OPTIONS = IDAQuadSetOptions(OLDOPTIONS,'NAME1',VALUE1,...) alters an 
%   existing options structure OLDOPTIONS.
%   
%   IDAQuadSetOptions with no input arguments displays all property names 
%   and their possible values.
%   
%IDAQuadSetOptions properties
%(See also the IDAS User Guide)
%
%ErrControl - Error control strategy for quadrature variables [ on | {off} ]
%   Specifies whether quadrature variables are included in the error test.
%RelTol - Relative tolerance for quadrature variables [ scalar {1e-4} ]
%   Specifies the relative tolerance for quadrature variables. This parameter is
%   used only if QuadErrCon=on.
%AbsTol - Absolute tolerance for quadrature variables [ scalar or vector {1e-6} ]
%   Specifies the absolute tolerance for quadrature variables. This parameter is
%   used only if QuadErrCon=on.
%
%SensDependent - Backward problem depending on sensitivities [ {false} | true ]
%   Specifies whether the backward problem quadrature right-hand side depends
%   on forward sensitivites. If TRUE, the right-hand side function provided for
%   this backward problem must have the appropriate type (see IDAQuadRhsFnB).
%
%   See also
%        IDAQuadInit, IDAQuadReInit.
%        IDAQuadInitB, IDAQuadReInitB

% Radu Serban <radu@llnl.gov>
% LLNS Start Copyright
% Copyright (c) 2013, Lawrence Livermore National Security
% This work was performed under the auspices of the U.S. Department 
% of Energy by Lawrence Livermore National Laboratory in part under 
% Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
% Produced at the Lawrence Livermore National Laboratory.
% All rights reserved.
% For details, see the LICENSE file.
% LLNS End Copyright
% $Revision: 1.1 $Date: 2006/08/10 18:01:03 $

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

options = idm_options(KeyNames,varargin{:});
