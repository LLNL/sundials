function options = CVodeSensSetOptions(varargin)
%CVodeSensSetOptions creates an options structure for FSA with CVODES.
%
%   Usage: OPTIONS = CVodeSensSetOptions('NAME1',VALUE1,'NAME2',VALUE2,...)
%          OPTIONS = CVodeSensSetOptions(OLDOPTIONS,'NAME1',VALUE1,...)
%
%   OPTIONS = CVodeSensSetOptions('NAME1',VALUE1,'NAME2',VALUE2,...) creates 
%   a CVODES options structure OPTIONS in which the named properties have 
%   the specified values. Any unspecified properties have default values. 
%   It is sufficient to type only the leading characters that uniquely 
%   identify the property. Case is ignored for property names. 
%   
%   OPTIONS = CVodeSensSetOptions(OLDOPTIONS,'NAME1',VALUE1,...) alters an 
%   existing options structure OLDOPTIONS.
%   
%   CVodeSensSetOptions with no input arguments displays all property names 
%   and their possible values.
%   
%CVodeSensSetOptions properties
%(See also the CVODES User Guide)
%
%method - FSA solution method [ 'Simultaneous' | {'Staggered'} ]
%   Specifies the FSA method for treating the nonlinear system solution for
%   sensitivity variables. In the simultaneous case, the nonlinear systems 
%   for states and all sensitivities are solved simultaneously. In the 
%   Staggered case, the nonlinear system for states is solved first and then
%   the nonlinear systems for all sensitivities are solved at the same time. 
%ParamField - Problem parameters  [ string ]
%   Specifies the name of the field in the user data structure (specified through
%   the 'UserData' field with CVodeSetOptions) in which the nominal values of the problem 
%   parameters are stored. This property is used only if  CVODES will use difference
%   quotient approximations to the sensitivity right-hand sides (see CVSensRhsFn).
%ParamList - Parameters with respect to which FSA is performed [ integer vector ]
%   Specifies a list of Ns parameters with respect to which sensitivities are to
%   be computed. This property is used only if CVODES will use difference-quotient
%   approximations to the sensitivity right-hand sides. Its length must be Ns, 
%   consistent with the number of columns of yS0 (see CVodeSensInit).
%ParamScales - Order of magnitude for problem parameters [ vector ]
%   Provides order of magnitude information for the parameters with respect to
%   which sensitivities are computed. This information is used if CVODES 
%   approximates the sensitivity right-hand sides or if CVODES estimates integration
%   tolerances for the sensitivity variables (see RelTol and AbsTol).
%RelTol - Relative tolerance for sensitivity variables [ positive scalar ]
%   Specifies the scalar relative tolerance for the sensitivity variables. 
%   See also AbsTol.
%AbsTol - Absolute tolerance for sensitivity variables [ row-vector or matrix ]
%   Specifies the absolute tolerance for sensitivity variables. AbsTol must be
%   either a row vector of dimension Ns, in which case each of its components is
%   used as a scalar absolute tolerance for the coresponding sensitivity vector,
%   or a N x Ns matrix, in which case each of its columns is used as a vector
%   of absolute tolerances for the corresponding sensitivity vector.
%   By default, CVODES estimates the integration tolerances for sensitivity 
%   variables, based on those for the states and on the order of magnitude 
%   information for the problem parameters specified through ParamScales.
%ErrControl - Error control strategy for sensitivity variables [ false | {true} ]
%   Specifies whether sensitivity variables are included in the error control test.
%   Note that sensitivity variables are always included in the nonlinear system
%   convergence test.
%DQtype - Type of DQ approx. of the sensi. RHS [{Centered} | Forward ]
%   Specifies whether to use centered (second-order) or forward (first-order)
%   difference quotient approximations of the sensitivity eqation right-hand 
%   sides. This property is used only if a user-defined sensitivity right-hand 
%   side function was not provided.
%DQparam - Cut-off parameter for the DQ approx. of the sensi. RHS [ scalar | {0.0} ]
%   Specifies the value which controls the selection of the difference-quotient 
%   scheme used in evaluating the sensitivity right-hand sides (switch between 
%   simultaneous or separate evaluations of the two components in the sensitivity 
%   right-hand side). The default value 0.0 indicates the use of simultaenous approximation
%   exclusively (centered or forward, depending on the value of DQtype.
%   For DQparam >= 1, CVODES uses a simultaneous approximation if the estimated
%   DQ perturbations for states and parameters are within a factor of DQparam, 
%   and separate approximations otherwise. Note that a value DQparam < 1
%   will inhibit switching! This property is used only if a user-defined sensitivity 
%   right-hand side function was not provided. 
%
%   See also
%        CVodeSensInit, CVodeSensReInit

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

if (nargin == 0) & (nargout == 0)
  fprintf('          method: [ Simultaneous | {Staggered} ]\n');
  fprintf('      ParamField: [ string ]\n');
  fprintf('       ParamList: [ integer vector ]\n');
  fprintf('     ParamScales: [ vector ]\n');
  fprintf('          RelTol: [ positive scalar ]\n');
  fprintf('          AbsTol: [ row-vector or matrix ]\n');
  fprintf('      ErrControl: [ false | {true} ]\n');
  fprintf('          DQtype: [ {Centered} | {Forward} ]\n');
  fprintf('         DQparam: [ scalar | {0.0} ]\n');
  fprintf('\n');
  return;
end

KeyNames = {
    'method'
    'ParamField'
    'ParamList'
    'ParamScales'
    'RelTol'
    'AbsTol'
    'ErrControl'
    'DQtype'
    'DQparam'
        };

options = cvm_options(KeyNames,varargin{:});

