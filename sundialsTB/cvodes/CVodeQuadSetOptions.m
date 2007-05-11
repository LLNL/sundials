function options = CVodeQuadSetOptions(varargin)
%CVodeQuadSetOptions creates an options structure for quadrature integration with CVODES.
%
%   Usage: OPTIONS = CVodeQuadSetOptions('NAME1',VALUE1,'NAME2',VALUE2,...)
%          OPTIONS = CVodeQuadSetOptions(OLDOPTIONS,'NAME1',VALUE1,...)
%          OPTIONS = CVodeQuadSetOptions(OLDOPTIONS,NEWOPTIONS)
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
%ErrControl - Error control strategy for quadrature variables [ on | {off} ]
%   Specifies whether quadrature variables are included in the error test.
%RelTol - Relative tolerance for quadrature variables [ scalar {1e-4} ]
%   Specifies the relative tolerance for quadrature variables. This parameter is
%   used only if ErrControl = on.
%AbsTol - Absolute tolerance for quadrature variables [ scalar or vector {1e-6} ]
%   Specifies the absolute tolerance for quadrature variables. This parameter is
%   used only if ErrControl = on.
%
%   See also
%        CVQuadInit, CVodeQuadReInit.
%        CVodeQuadInitB, CVodeQuadReInitB

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/08/10 17:59:57 $

% If called without input and output arguments, print out the possible keywords

if (nargin == 0) && (nargout == 0)
  fprintf('      ErrControl: [ on | {off} ]\n');
  fprintf('          RelTol: [ positive scalar {1e-4} ]\n');
  fprintf('          AbsTol: [ positive scalar or vector {1e-6} ]\n');
  fprintf('\n');
  return;
end

KeyNames = {
    'ErrControl'
    'RelTol'
    'AbsTol'
    };

options = cvm_options(KeyNames,varargin{:});
