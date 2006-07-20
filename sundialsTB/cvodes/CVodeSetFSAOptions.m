function options = CVodeSetFSAOptions(varargin)
%CVodeSetFSAOptions creates an options structure for FSA with CVODES.
%
%   Usage: OPTIONS = CVodeSetFSAOptions('NAME1',VALUE1,'NAME2',VALUE2,...)
%          OPTIONS = CVodeSetFSAOptions(OLDOPTIONS,'NAME1',VALUE1,...)
%          OPTIONS = CVodeSetFSAOptions(OLDOPTIONS,NEWOPTIONS)
%
%   OPTIONS = CVodeSetFSAOptions('NAME1',VALUE1,'NAME2',VALUE2,...) creates 
%   a CVODES options structure OPTIONS in which the named properties have 
%   the specified values. Any unspecified properties have default values. 
%   It is sufficient to type only the leading characters that uniquely 
%   identify the property. Case is ignored for property names. 
%   
%   OPTIONS = CVodeSetFSAOptions(OLDOPTIONS,'NAME1',VALUE1,...) alters an 
%   existing options structure OLDOPTIONS.
%   
%   OPTIONS = CVodeSetFSAOptions(OLDOPTIONS,NEWOPTIONS) combines an existing 
%   options structure OLDOPTIONS with a new options structure NEWOPTIONS. 
%   Any new properties overwrite corresponding old properties. 
%   
%   CVodeSetFSAOptions with no input arguments displays all property names 
%   and their possible values.
%   
%CVodeSetFSAOptions properties
%(See also the CVODES User Guide)
%
%ParamField - Problem parameters  [ string ]
%   Specifies the name of the field in the user data structure (passed as an 
%   argument to CVodeMalloc) in which the nominal values of the problem 
%   parameters are stored. This property is used only if  CVODES will use difference
%   quotient approximations to the sensitivity right-hand sides (see SensRhsFn).
%ParamList - Parameters with respect to which FSA is performed [ integer vector ]
%   Specifies a list of Ns parameters with respect to which sensitivities are to
%   be computed. This property is used only if CVODES will use difference-quotient
%   approximations to the sensitivity right-hand sides (see SensRhsFn below).
%   Its length must be Ns, consistent with the number of columns of yS0
%   (see CVodeSensMalloc).
%ParamScales - Order of magnitude for problem parameters [ vector ]
%   Provides order of magnitude information for the parameters with respect to
%   which sensitivities are computed. This information is used if CVODES 
%   approximates the sensitivity right-hand sides (see SensRhsFn below) or if CVODES 
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
%   By default, CVODES estimates the integration tolerances for sensitivity 
%   variables, based on those for the states and on the order of magnitude 
%   information for the problem parameters specified through ParamScales.
%SensErrControl - Error control strategy for sensitivity variables [ on | {off} ]
%   Specifies whether sensitivity variables are included in the error control test.
%   Note that sensitivity variables are always included in the nonlinear system
%   convergence test.
%SensRhsFn - Sensitivity right-hand side function [ function ]
%   Specifies a user-supplied function to evaluate the sensitivity right-hand 
%   sides. If not specified, CVODES uses a default internal difference-quotient 
%   function to approximate the sensitivity right-hand sides.
%SensDQtype - Type of DQ approx. of the sensi. RHS [{Centered} | Forward ]
%   Specifies whether to use centered (second-order) or forward (first-order)
%   difference quotient approximations of the sensitivity eqation right-hand 
%   sides. This property is used only if a user-defined sensitivity right-hand 
%   side function was not provided.
%SensDQparam - Cut-off parameter for the DQ approx. of the sensi. RHS [ scalar | {0.0} ]
%   Specifies the value which controls the selection of the difference-quotient 
%   scheme used in evaluating the sensitivity right-hand sides (switch between 
%   simultaneous or separate evaluations of the two components in the sensitivity 
%   right-hand side). The default value 0.0 indicates the use of simultaenous approximation
%   exclusively (centered or forward, depending on the value of SensDQtype.
%   For SensDQparam >= 1, CVODES uses a simultaneous approximation if the estimated
%   DQ perturbations for states and parameters are within a factor of SensDQparam, 
%   and separate approximations otherwise. Note that a value SensDQparam < 1
%   will inhibit switching! This property is used only if a user-defined sensitivity 
%   right-hand side function was not provided. 
%
%   See also
%        CVodeSensMalloc, CVSensRhsFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/03/07 01:19:50 $

% Based on Matlab's ODESET function

% Print out possible values of properties.
if (nargin == 0) & (nargout == 0)
  fprintf('      ParamField: [ string ]\n');
  fprintf('       ParamList: [ integer vector ]\n');
  fprintf('     ParamScales: [ vector ]\n');
  fprintf('      SensRelTol: [ positive scalar ]\n');
  fprintf('      SensAbsTol: [ row-vector or matrix ]\n');
  fprintf('  SensErrControl: [ off | {on} ]\n');
  fprintf('       SensRhsFn: [ function ]\n');
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
    'SensRhsFn       '
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
                     'or an options structure\ncreated with CVodeSetFSAOptions.'], i));
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
