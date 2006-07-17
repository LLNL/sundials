function options = IDASetOptions(varargin)
%IDASetOptions creates an options structure for IDAS.
%
%   Usage: OPTIONS = IDASetOptions('NAME1',VALUE1,'NAME2',VALUE2,...)
%          OPTIONS = IDASetOptions(OLDOPTIONS,'NAME1',VALUE1,...)
%          OPTIONS = IDASetOptions(OLDOPTIONS,NEWOPTIONS)
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
%   OPTIONS = IDASetOptions(OLDOPTIONS,NEWOPTIONS) combines an existing 
%   options structure OLDOPTIONS with a new options structure NEWOPTIONS. 
%   Any new properties overwrite corresponding old properties. 
%   
%   IDASetOptions with no input arguments displays all property names 
%   and their possible values.
%   
%IDASetOptions properties
%(See also the IDAS User Guide)
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
%VariableTypes - Alg./diff. variables [ vector ]
%ConstraintTypes - Simple bound constraints [ vector ]
%
%ICcalculation - Consistent IC calculation [{None}|FindAlgebraic|FindAll]
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
%Quadratures - Quadrature integration [ on | {off} ]
%   Enables or disables quadrature integration.
%QuadRhsFn - Quadrature residual function [ function ]
%   Specifies the user-supplied function to evaluate the integrand for
%   quadrature computations. See IDAQuadRhsfn.
%QuadInitCond - Initial conditions for quadrature variables [ vector ]
%   Specifies the initial conditions for quadrature variables.
%QuadErrControl - Error control strategy for quadrature variables [ on | {off} ]
%   Specifies whether quadrature variables are included in the error test.
%QuadRelTol - Relative tolerance for quadrature variables [ scalar {1e-4} ]
%   Specifies the relative tolerance for quadrature variables. This parameter is
%   used only if QuadErrCon=on.
%QuadAbsTol - Absolute tolerance for quadrature variables [ scalar or vector {1e-6} ]
%   Specifies the absolute tolerance for quadrature variables. This parameter is
%   used only if QuadErrCon=on.
%
%ASANumDataPoints - Number of data points for ASA [ integer | {100} ]
%   Specifies the (maximum) number of integration steps between two consecutive
%   check points.
%ASAInterpType - Type of interpolation [ Polynomial | {Hermite} ]
%   Specifies the type of interpolation used for estimating the forward solution
%   during the backward integration phase. At this time, the only option is
%   'Hermite', specifying cubic Hermite interpolation.
%
%MonitorFn - User-provied monitoring function [ function ]
%   Specifies a function that is called after each successful integration step.
%   This function must have type IDAMonitorFn. A simple monitoring function,
%   IDAMonitor is provided with IDAS.
%MonitorData - User-provied data for the monitoring function [ struct ]
%   Specifies a data structure that is passed to the Monitor function every time
%   it is called. 
%
%   See also
%        IDARootFn, IDAQuadRhsFn
%        IDADenseJacFn, IDABandJacFn, IDAJacTimesVecFn
%        IDAPrecSetupFn, IDAPrecSolveFn
%        IDAGlocalFn, IDAGcommFn
%        IDAMonitorFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/03/15 19:31:25 $

% Based on Matlab's ODESET function

% Print out possible values of properties.
if (nargin == 0) & (nargout == 0)
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
  fprintf('   VariableTypes: [ vector ]\n');
  fprintf(' ConstraintTypes: [ vector ]\n');
  fprintf('\n');
  fprintf('   ICcalculation: [ {None} | FindAlgebraic | FindAll ]\n');
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
  fprintf('     Quadratures: [ on | {off} ]\n');
  fprintf('       QuadRhsFn: [ function ]\n');
  fprintf('    QuadInitCond: [ vector ]\n');
  fprintf('  QuadErrControl: [ on | {off} ]\n');
  fprintf('      QuadRelTol: [ positive scalar {1e-4} ]\n');
  fprintf('      QuadAbsTol: [ positive scalar or vector {1e-6} ]\n');
  fprintf('\n');
  fprintf('    ASANumPoints: [ integer | {100} ]\n');
  fprintf('   ASAInterpType: [ Polynomial | {Hermite} ]\n');
  fprintf('\n');
  fprintf('       MonitorFn: [ function ]\n');
  fprintf('     MonitorData: [ struct ]\n');
  fprintf('\n');
  return;
end

Names = [
    'RelTol          '
    'AbsTol          '
    'MaxNumSteps     '
    'InitialStep     '
    'MaxStep         '
    'MaxOrder        '
    'StopTime        '
    'RootsFn         '
    'NumRoots        '
    'VariableTypes   '
    'ConstraintTypes '
    'ICcalculation   '
    'LinearSolver    '
    'JacobianFn      '
    'PrecModule      '
    'PrecSetupFn     '
    'PrecSolveFn     '
    'KrylovMaxDim    '
    'GramSchmidtType '
    'GlocalFn        '
    'GcommFn         '
    'LowerBwidth     '
    'UpperBwidth     '
    'LowerBwidthDQ   '
    'UpperBwidthDQ   '
    'Quadratures     '
    'QuadRhsFn       '
    'QuadInitCond    '
    'QuadErrControl  '
    'QuadRelTol      '
    'QuadAbsTol      '
    'ASANumPoints    '
    'ASAInterpType   '
    'MonitorFn       '
    'MonitorData     '
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
                     'or an options structure\ncreated with IDASetOptions.'], i));
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
