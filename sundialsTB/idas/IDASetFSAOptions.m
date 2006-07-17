function options = IDASetFSAOptions(varargin)
%IDASetFSAOptions creates an options structure for FSA with IDAS.
%
%   Usage: OPTIONS = IDASetFSAOptions('NAME1',VALUE1,'NAME2',VALUE2,...)
%          OPTIONS = IDASetFSAOptions(OLDOPTIONS,'NAME1',VALUE1,...)
%          OPTIONS = IDASetFSAOptions(OLDOPTIONS,NEWOPTIONS)
%
%   OPTIONS = IDASetFSAOptions('NAME1',VALUE1,'NAME2',VALUE2,...) creates 
%   a IDAS options structure OPTIONS in which the named properties have 
%   the specified values. Any unspecified properties have default values. 
%   It is sufficient to type only the leading characters that uniquely 
%   identify the property. Case is ignored for property names. 
%   
%   OPTIONS = IDASetFSAOptions(OLDOPTIONS,'NAME1',VALUE1,...) alters an 
%   existing options structure OLDOPTIONS.
%   
%   OPTIONS = IDASetFSAOptions(OLDOPTIONS,NEWOPTIONS) combines an existing 
%   options structure OLDOPTIONS with a new options structure NEWOPTIONS. 
%   Any new properties overwrite corresponding old properties. 
%   
%   IDASetFSAOptions with no input arguments displays all property names 
%   and their possible values.
%   
%IDASetFSAOptions properties
%(See also the IDAS User Guide)
%
%FSAParamField - Problem parameters  [ string ]
%   Specifies the name of the field in the user data structure (passed as an 
%   argument to IDAMalloc) in which the nominal values of the problem 
%   parameters are stored. This property is used only if  IDAS will use difference
%   quotient approximations to the sensitivity residuals (see SensRHS and 
%   SensRHStype).
%FSAParamList - Parameters with respect to which FSA is performed [ integer vector ]
%   Specifies a list of Ns parameters with respect to which sensitivities are to
%   be computed. This property is used only if IDAS will use difference-quotient
%   approximations to the sensitivity residuals (see SensRHS and SensRHStype). 
%   Its length must be Ns, consistent with the number of columns of FSAinitCond.
%FSAParamScales - Order of magnitude for problem parameters [ vector ]
%   Provides order of magnitude information for the parameters with respect to
%   which sensitivities are computed. This information is used if IDAS 
%   approximates the sensitivity residuals (see SensRHS) or if IDAS estimates 
%   integration tolerances for the sensitivity variables (see FSAReltol and FSAAbsTol).
%FSARelTol - Relative tolerance for sensitivity variables [ positive scalar ]
%   Specifies the scalar relative tolerance for the sensitivity variables. 
%   See FSAAbsTol.
%FSAAbsTol - Absolute tolerance for sensitivity variables [ row-vector or matrix ]
%   Specifies the absolute tolerance for sensitivity variables. FSAAbsTol must be
%   either a row vector of dimension Ns, in which case each of its components is
%   used as a scalar absolute tolerance for the coresponding sensitivity vector,
%   or a N x Ns matrix, in which case each of its columns is used as a vector
%   of absolute tolerances for the corresponding sensitivity vector.
%   By default, IDAS estimates the integration tolerances for sensitivity 
%   variables, based on those for the states and on the order of magnitude 
%   information for the problem parameters specified through ParamScales.
%FSAErrControl - Error control strategy for sensitivity variables [ on | {off} ]
%   Specifies whether sensitivity variables are included in the error control test.
%   Note that sensitivity variables are always included in the nonlinear system
%   convergence test.
%FSAResFn - Sensitivity residual function [ function ]
%   Specifies a user-supplied function to evaluate the sensitivity residuals.
%   See SensRHStype. By default, IDAS uses an internal difference-quotient 
%   function to approximate the sensitivity residuals.
%FSADQparam - Parameter for the DQ approx. of the sensi. residual [ scalar | {0.0} ]
%   Specifies the value which controls the selection of the difference-quotient 
%   scheme used in evaluating the sensitivity residuals. This property is 
%   used only if IDAS will use difference-quotient approximations. The default 
%   value 0.0 indicates the use of the second-order centered directional derivative 
%   formula exclusively. Otherwise, the magnitude of FSADQparam and its sign 
%   (positive or negative) indicates whether this switching is done with regard 
%   to (centered or forward) finite differences, respectively.
%
%   See also
%        IDASensRhsFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/03/07 01:19:50 $

% Based on Matlab's ODESET function

% Print out possible values of properties.
if (nargin == 0) & (nargout == 0)
  fprintf('   FSAParamField: [ string ]\n');
  fprintf('    FSAParamList: [ integer vector ]\n');
  fprintf('  FSAParamScales: [ vector ]\n');
  fprintf('       FSARelTol: [ positive scalar ]\n');
  fprintf('       FSAAbsTol: [ row-vector or matrix ]\n');
  fprintf('   FSAErrControl: [ off | {on} ]\n');
  fprintf('        FSAResFn: [ function ]\n');
  fprintf('      FSADQparam: [ scalar | {0.0} ]\n');
  return;
end

Names = [
    'FSAParamField   '
    'FSAParamList    '
    'FSAParamScales  '
    'FSARelTol       '
    'FSAAbsTol       '
    'FSAErrControl   '
    'FSAResFn        '
    'FSADQparam      '
        ];
[m,n] = size(Names);
names = lower(Names);

% Combine all leading options structures o1, o2, ... in (o1,o2,...).
options = [];
for j = 1:m
  options.(deblank(Names(j,:))) = [];
end
i = 1;
while i <= nargin
  arg = varargin{i};
  if isstr(arg)                         % arg is an option name
    break;
  end
  if ~isempty(arg)                      % [] is a valid options argument
    if ~isa(arg,'struct')
      error(sprintf(['Expected argument %d to be a string property name ' ...
                     'or an options structure\ncreated with IDASetFSAOptions.'], i));
    end
    for j = 1:m
      if any(strcmp(fieldnames(arg),deblank(Names(j,:))))
        val = arg.(deblank(Names(j,:)));
      else
        val = [];
      end
      if ~isempty(val)
        options.(deblank(Names(j,:))) = val;
      end
    end
  end
  i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
  error('Arguments must occur in name-value pairs.');
end
expectval = 0;                          % start expecting a name, not a value
while i <= nargin
  arg = varargin{i};
    
  if ~expectval
    if ~isstr(arg)
      error(sprintf('Expected argument %d to be a string property name.', i));
    end
    
    lowArg = lower(arg);
    j = strmatch(lowArg,names);
    if isempty(j)                       % if no matches
      error(sprintf('Unrecognized property name ''%s''.', arg));
    elseif length(j) > 1                % if more than one match
      % Check for any exact matches (in case any names are subsets of others)
      k = strmatch(lowArg,names,'exact');
      if length(k) == 1
        j = k;
      else
        msg = sprintf('Ambiguous property name ''%s'' ', arg);
        msg = [msg '(' deblank(Names(j(1),:))];
        for k = j(2:length(j))'
          msg = [msg ', ' deblank(Names(k,:))];
        end
        msg = sprintf('%s).', msg);
        error(msg);
      end
    end
    expectval = 1;                      % we expect a value next
    
  else
    options.(deblank(Names(j,:))) = arg;
    expectval = 0;
      
  end
  i = i + 1;
end

if expectval
  error(sprintf('Expected value for property ''%s''.', arg));
end
