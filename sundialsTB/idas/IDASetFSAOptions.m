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
%ParamField - Problem parameters  [ string ]
%   Specifies the name of the field in the user data structure (passed as an 
%   argument to IDAMalloc) in which the nominal values of the problem 
%   parameters are stored. This property is used only if  IDAS will use difference
%   quotient approximations to the sensitivity residuals (see SensResFn).
%ParamList - Parameters with respect to which FSA is performed [ integer vector ]
%   Specifies a list of Ns parameters with respect to which sensitivities are to
%   be computed. This property is used only if IDAS will use difference-quotient
%   approximations to the sensitivity residuals (see SensResFn below).
%   Its length must be Ns, consistent with the number of columns of yS0
%   (see IDASensMalloc).
%ParamScales - Order of magnitude for problem parameters [ vector ]
%   Provides order of magnitude information for the parameters with respect to
%   which sensitivities are computed. This information is used if IDAS 
%   approximates the sensitivity residuals (see SensResFn below) or if IDAS 
%   estimates integration tolerances for the sensitivity variables (see SensReltol 
%   and SensAbsTol).
%SensRelTol - Relative tolerance for sensitivity variables [ positive scalar ]
%   Specifies the scalar relative tolerance for the sensitivity variables. 
%   See also SensAbsTol.
%SensAbsTol - Absolute tolerance for sensitivity variables [ row-vector or matrix ]
%   Specifies the absolute tolerance for sensitivity variables. SensAbsTol must be
%   either a row vector of dimension Ns, in which case each of its components is
%   used as a scalar absolute tolerance for the coresponding sensitivity vector,
%   or a N x Ns matrix, in which case each of its columns is used as a vector
%   of absolute tolerances for the corresponding sensitivity vector.
%   By default, IDAS estimates the integration tolerances for sensitivity 
%   variables, based on those for the states and on the order of magnitude 
%   information for the problem parameters specified through ParamScales.
%SensErrControl - Error control strategy for sensitivity variables [ on | {off} ]
%   Specifies whether sensitivity variables are included in the error control test.
%   Note that sensitivity variables are always included in the nonlinear system
%   convergence test.
%SensResFn - Sensitivity residual function [ function ]
%   Specifies a user-supplied function to evaluate the sensitivity right-hand 
%   sides. If not specified, IDAS uses a default internal difference-quotient 
%   function to approximate the sensitivity residuals.
%SensDQtype - Type of DQ approx. of the sensi. RHS [{Centered} | Forward ]
%   Specifies whether to use centered (second-order) or forward (first-order)
%   difference quotient approximations of the sensitivity eqation right-hand 
%   sides. This property is used only if a user-defined sensitivity right-hand 
%   side function was not provided.
%SensDQparam - Cut-off parameter for the DQ approx. of the sensi. RHS [ scalar | {0.0} ]
%   Specifies the value which controls the selection of the difference-quotient 
%   scheme used in evaluating the sensitivity residuals (switch between 
%   simultaneous or separate evaluations of the two components in the sensitivity 
%   residual). The default value 0.0 indicates the use of simultaenous approximation
%   exclusively (centered or forward, depending on the value of SensDQtype.
%   For SensDQparam >= 1, IDAS uses a simultaneous approximation if the estimated
%   DQ perturbations for states and parameters are within a factor of SensDQparam, 
%   and separate approximations otherwise. Note that a value SensDQparam < 1
%   will inhibit switching! This property is used only if a user-defined sensitivity 
%   residual function was not provided. 
%
%   See also
%        IDASensMalloc, IDASensResFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/07/17 16:49:50 $

% Based on Matlab's ODESET function

% Print out possible values of properties.
if (nargin == 0) & (nargout == 0)
  fprintf('      ParamField: [ string ]\n');
  fprintf('       ParamList: [ integer vector ]\n');
  fprintf('     ParamScales: [ vector ]\n');
  fprintf('      SensRelTol: [ positive scalar ]\n');
  fprintf('      SensAbsTol: [ row-vector or matrix ]\n');
  fprintf('  SensErrControl: [ off | {on} ]\n');
  fprintf('       SensResFn: [ function ]\n');
  fprintf('      SensDQtype: [ {Centered} | {Forward} ]\n');
  fprintf('     SensDQparam: [ scalar | {0.0} ]\n');
  return;
end

Names = [
    'ParamField      '
    'ParamList       '
    'ParamScales     '
    'SensRelTol      '
    'SensAbsTol      '
    'SensErrControl  '
    'SensResFn       '
    'SensDQtype      '
    'SensDQparam     '
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
