function options = IDASetOptions(varargin)
%IDASetOptions creates an options structure for IDAS.
%
%   Usage: OPTIONS = IDASetOptions('NAME1',VALUE1,'NAME2',VALUE2,...)
%          OPTIONS = IDASetOptions(OLDOPTIONS,'NAME1',VALUE1,...)
%
%   OPTIONS = IDASetOptions('NAME1',VALUE1,'NAME2',VALUE2,...) creates 
%   a IDAS options structure OPTIONS in which the named properties have 
%   the specified values. Any unspecified properties have default values. 
%   It is sufficient to type only the leading characters that uniquely 
%   identify the property. Case is ignored for property names. 
%   
%   OPTIONS = IDASetOptions(OLDOPTIONS,'NAME1',VALUE1,...) alters an 
%   existing options structure OLDOPTIONS.
%   
%   IDASetOptions with no input arguments displays all property names 
%   and their possible values.
%   
%IDASetOptions properties
%(See also the IDAS User Guide)
%
%UserData - User data passed unmodified to all functions [ empty ]
%   If UserData is not empty, all user provided functions will be
%   passed the problem data as their last input argument. For example,
%   the RES function must be defined as R = DAEFUN(T,YY,TP,DATA).
%
%RelTol - Relative tolerance [ positive scalar | {1e-4} ]
%   RelTol defaults to 1e-4 and is applied to all components of the solution
%   vector. See AbsTol.
%AbsTol - Absolute tolerance [ positive scalar or vector | {1e-6} ]
%   The relative and absolute tolerances define a vector of error weights
%   with components
%     ewt(i) = 1/(RelTol*|y(i)| + AbsTol)    if AbsTol is a scalar
%     ewt(i) = 1/(RelTol*|y(i)| + AbsTol(i)) if AbsTol is a vector
%   This vector is used in all error and convergence tests, which
%   use a weighted RMS norm on all error-like vectors v:
%     WRMSnorm(v) = sqrt( (1/N) sum(i=1..N) (v(i)*ewt(i))^2 ),
%   where N is the problem dimension.
%MaxNumSteps - Maximum number of steps [positive integer | {500}]
%   IDASolve will return with an error after taking MaxNumSteps internal steps
%   in its attempt to reach the next output time.
%InitialStep - Suggested initial stepsize [ positive scalar ]
%   By default, IDASolve estimates an initial stepsize h0 at the initial time 
%   t0 as the solution of
%     WRMSnorm(h0^2 ydd / 2) = 1
%   where ydd is an estimated second derivative of y(t0).
%MaxStep - Maximum stepsize [ positive scalar | {inf} ]
%   Defines an upper bound on the integration step size.
%MaxOrder - Maximum method order [ 1-5 for BDF | {5} ]
%   Defines an upper bound on the linear multistep method order.
%StopTime - Stopping time [ scalar ]
%   Defines a value for the independent variable past which the solution 
%   is not to proceed. 
%RootsFn - Rootfinding function [ function ]
%   To detect events (roots of functions), set this property to the event 
%   function. See IDARootFn.
%NumRoots - Number of root functions [ integer | {0} ]
%   Set NumRoots to the number of functions for which roots are monitored.
%   If NumRoots is 0, rootfinding is disabled.
%
%SuppressAlgVars - Suppres algebraic vars. from error test [ on | {off} ]
%VariableTypes - Alg./diff. variables [ vector ]
%ConstraintTypes - Simple bound constraints [ vector ]
%
%LinearSolver - Linear solver type [{Dense}|Band|GMRES|BiCGStab|TFQMR]
%   Specifies the type of linear solver to be used for the Newton nonlinear 
%   solver. Valid choices are: Dense (direct, dense Jacobian), Band (direct, 
%   banded Jacobian), GMRES (iterative, scaled preconditioned GMRES), 
%   BiCGStab (iterative, scaled preconditioned stabilized BiCG), TFQMR 
%   (iterative, scaled transpose-free QMR).
%   The GMRES, BiCGStab, and TFQMR are matrix-free linear solvers.
%JacobianFn - Jacobian function [ function ]
%   This propeerty is overloaded. Set this value to a function that returns 
%   Jacobian information consistent with the linear solver used (see Linsolver). 
%   If not specified, IDAS uses difference quotient approximations. 
%   For the Dense linear solver, JacobianFn must be of type IDADenseJacFn and 
%   must return a dense Jacobian matrix. For the Band linear solver, JacobianFn 
%   must be of type IDABandJacFn and must return a banded Jacobian matrix. 
%   For the iterative linear solvers, GMRES, BiCGStab, and TFQMR, JacobianFn must 
%   be of type IDAJacTimesVecFn and must return a Jacobian-vector product.
%KrylovMaxDim - Maximum number of Krylov subspace vectors [ integer | {5} ]
%   Specifies the maximum number of vectors in the Krylov subspace. This property 
%   is used only if an iterative linear solver, GMRES, BiCGStab, or TFQMR is used 
%   (see LinSolver).
%GramSchmidtType - Gram-Schmidt orthogonalization [ Classical | {Modified} ]
%   Specifies the type of Gram-Schmidt orthogonalization (classical or modified).
%   This property is used only if the GMRES linear solver is used (see LinSolver).
%PrecModule - Preconditioner module [ BBDPre | {UserDefined} ]
%   If PrecModule = 'UserDefined', then the user must provide at least a 
%   preconditioner solve function (see PrecSolveFn)
%   IDAS provides one general-purpose preconditioner module, BBDPre, which can 
%   be only used with parallel vectors. It provide a preconditioner matrix that 
%   is block-diagonal with banded blocks. The blocking corresponds to the 
%   distribution of the dependent variable vector y among the processors. 
%   Each preconditioner block is generated from the Jacobian of the local part 
%   (on the current processor) of a given function g(t,y,yp) approximating
%   f(t,y,yp) (see GlocalFn). The blocks are generated by a difference quotient 
%   scheme on each processor independently. This scheme utilizes an assumed 
%   banded structure with given half-bandwidths, mldq and mudq (specified through 
%   LowerBwidthDQ and UpperBwidthDQ, respectively). However, the banded Jacobian 
%   block kept by the scheme has half-bandwiths ml and mu (specified through 
%   LowerBwidth and UpperBwidth), which may be smaller.
%PrecSetupFn - Preconditioner setup function [ function ]
%   If PrecType is not 'None', PrecSetupFn specifies an optional function which,
%   together with PrecSolve, defines the preconditioner matrix, which must be an
%   aproximation to the Newton matrix. PrecSetupFn must be of type IDAPrecSetupFn. 
%PrecSolveFn - Preconditioner solve function [ function ]
%   If PrecType is not 'None', PrecSolveFn specifies a required function which 
%   must solve a linear system Pz = r, for given r. PrecSolveFn must be of type
%   IDAPrecSolveFn.
%GlocalFn - Local residual approximation function for BBDPre [ function ]
%   If PrecModule is BBDPre, GlocalFn specifies a required function that
%   evaluates a local approximation to the DAE residual. GlocalFn must
%   be of type IDAGlocFn.
%GcommFn - Inter-process communication function for BBDPre [ function ]
%   If PrecModule is BBDPre, GcommFn specifies an optional function
%   to perform any inter-process communication required for the evaluation of
%   GlocalFn. GcommFn must be of type IDAGcommFn.
%LowerBwidth - Jacobian/preconditioner lower bandwidth [ integer | {0} ]
%   This property is overloaded. If the Band linear solver is used (see LinSolver),
%   it specifies the lower half-bandwidth of the band Jacobian approximation. 
%   If one of the three iterative linear solvers, GMRES, BiCGStab, or TFQMR is used 
%   (see LinSolver) and if the BBDPre preconditioner module in IDAS is used 
%   (see PrecModule), it specifies the lower half-bandwidth of the retained 
%   banded approximation of the local Jacobian block.
%   LowerBwidth defaults to 0 (no sub-diagonals).
%UpperBwidth - Jacobian/preconditioner upper bandwidth [ integer | {0} ]
%   This property is overloaded. If the Band linear solver is used (see LinSolver),
%   it specifies the upper half-bandwidth of the band Jacobian approximation. 
%   If one of the three iterative linear solvers, GMRES, BiCGStab, or TFQMR is used 
%   (see LinSolver) and if the BBDPre preconditioner module in IDAS is used 
%   (see PrecModule), it specifies the upper half-bandwidth of the retained 
%   banded approximation of the local Jacobian block.
%   UpperBwidth defaults to 0 (no super-diagonals).
%LowerBwidthDQ - BBDPre preconditioner DQ lower bandwidth [ integer | {0} ]
%   Specifies the lower half-bandwidth used in the difference-quotient Jacobian
%   approximation for the BBDPre preconditioner (see PrecModule).
%UpperBwidthDQ - BBDPre preconditioner DQ upper bandwidth [ integer | {0} ]
%   Specifies the upper half-bandwidth used in the difference-quotient Jacobian
%   approximation for the BBDPre preconditioner (see PrecModule).
%
%MonitorFn - User-provied monitoring function [ function ]
%   Specifies a function that is called after each successful integration step.
%   This function must have type IDAMonitorFn or IDAMonitorFnB, depending on
%   whether these options are for a forward or a backward problem, respectively.
%   Sample monitoring functions IDAMonitor and IDAMonitorB are provided
%   with IDAS.
%MonitorData - User-provied data for the monitoring function [ struct ]
%   Specifies a data structure that is passed to the MonitorFn function every time
%   it is called. 
%
%SensDependent - Backward problem depending on sensitivities [ {false} | true ]
%   Specifies whether the backward problem right-hand side depends on 
%   forward sensitivites. If SUNTRUE, the residual function provided for
%   this backward problem must have the appropriate type (see IDAResFnB).
%
%ErrorMessages - Post error/warning messages [ {true} | false ]
%   Note that any errors in IDAInit will result in a Matlab error, thus
%   stoping execution. Only subsequent calls to IDAS functions will respect
%   the value specified for 'ErrorMessages'.
%
%NOTES:
%
%   The properties listed above that can only be used for forward problems
%   are: ConstraintTypes, StopTime, RootsFn, and NumRoots.
%
%   The property SensDependent is relevant only for backward problems.
%
%   See also
%        IDAInit, IDAReInit, IDAInitB, IDAReInitB
%        IDAResFn, IDARootFn
%        IDADenseJacFn, IDABandJacFn, IDAJacTimesVecFn
%        IDAPrecSetupFn, IDAPrecSolveFn
%        IDAGlocalFn, IDAGcommFn
%        IDAMonitorFn
%        IDAResFnB
%        IDADenseJacFnB, IDABandJacFnB, IDAJacTimesVecFnB
%        IDAPrecSetupFnB, IDAPrecSolveFnB
%        IDAGlocalFnB, IDAGcommFnB
%        IDAMonitorFnB

% Radu Serban <radu@llnl.gov>
% SUNDIALS Copyright Start
% Copyright (c) 2002-2021, Lawrence Livermore National Security
% and Southern Methodist University.
% All rights reserved.
%
% See the top-level LICENSE and NOTICE files for details.
%
% SPDX-License-Identifier: BSD-3-Clause
% SUNDIALS Copyright End
% $Revision$Date: 2007/08/21 17:38:43 $


% If called without input and output arguments, print out the possible keywords

if (nargin == 0) & (nargout == 0)
  fprintf('        UserData: [ empty ]\n');
  fprintf('\n');
  fprintf('          RelTol: [ positive scalar | {1e-4} ]\n');
  fprintf('          AbsTol: [ positive scalar or vector | {1e-6} ]\n');
  fprintf('     MaxNumSteps: [ positive integer | {500} ]\n');
  fprintf('     InitialStep: [ positive scalar ]\n');
  fprintf('         MaxStep: [ positive scalar | {inf} ]\n');
  fprintf('        MaxOrder: [ 1-12 for Adams, 1-5 for BDF | {5} ]\n');
  fprintf('        StopTime: [ scalar ]\n');
  fprintf('         RootsFn: [ function ]\n');
  fprintf('        NumRoots: [ integer | {0} ]\n');
  fprintf('\n');
  fprintf(' SuppressAlgVars: [ on | {off} ]\n');
  fprintf('   VariableTypes: [ vector ]\n');
  fprintf(' ConstraintTypes: [ vector ]\n');
  fprintf('\n');
  fprintf('    LinearSolver: [ {Dense} | Band | GMRES | BiCGStab | TFQMR ]\n');
  fprintf('      JacobianFn: [ function ]\n');
  fprintf('    KrylovMaxDim: [ integer | {5} ]\n');
  fprintf(' GramSchmidtType: [ Classical | {Modified} ]\n');
  fprintf('      PrecModule: [ BBDPre | {UserDefined} ]\n');
  fprintf('     PrecSetupFn: [ function ]\n');
  fprintf('     PrecSolveFn: [ function ]\n');
  fprintf('        GlocalFn: [ function ]\n');
  fprintf('         GcommFn: [ function ]\n');
  fprintf('     LowerBwidth: [ integer | {0} ]\n');
  fprintf('     UpperBwidth: [ integer | {0} ]\n');
  fprintf('   LowerBwidthDQ: [ integer | {0} ]\n');
  fprintf('   UpperBwidthDQ: [ integer | {0} ]\n');
  fprintf('\n');
  fprintf('       MonitorFn: [ function ]\n');
  fprintf('     MonitorData: [ struct ]\n');
  fprintf('\n');
  fprintf('   SensDependent: [ {false} | true ]\n');
  fprintf('\n');
  fprintf('   ErrorMessages: [ false | {true} ]\n');
  fprintf('\n');
  return;
end

KeyNames = {
    'UserData'
    'RelTol'
    'AbsTol'
    'MaxNumSteps'
    'InitialStep'
    'MaxStep'
    'MaxOrder'
    'StopTime'
    'RootsFn'
    'NumRoots'
    'VariableTypes'
    'ConstraintTypes'
    'SuppressAlgVars'
    'LinearSolver'
    'JacobianFn'
    'PrecModule'
    'PrecSetupFn'
    'PrecSolveFn'
    'KrylovMaxDim'
    'GramSchmidtType'
    'GlocalFn'
    'GcommFn'
    'LowerBwidth'
    'UpperBwidth'
    'LowerBwidthDQ'
    'UpperBwidthDQ'
    'MonitorFn'
    'MonitorData'
    'SensDependent'
    'ErrorMessages'
           };

options = idm_options(KeyNames,varargin{:});