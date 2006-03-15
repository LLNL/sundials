function options = KINSetOptions(varargin)
%KINSetOptions creates an options structure for KINSOL.
%
%   Usage:
%
%   options = KINSetOptions('NAME1',VALUE1,'NAME2',VALUE2,...) creates a KINSOL
%   options structure options in which the named properties have the
%   specified values. Any unspecified properties have default values. It is
%   sufficient to type only the leading characters that uniquely identify the
%   property. Case is ignored for property names. 
%   
%   options = KINSetOptions(oldoptions,'NAME1',VALUE1,...) alters an existing 
%   options structure oldoptions.
%   
%   options = KINSetOptions(oldoptions,newoptions) combines an existing options 
%   structure oldoptions with a new options structure newoptions. Any new 
%   properties overwrite corresponding old properties. 
%   
%   KINSetOptions with no input arguments displays all property names and their
%   possible values.
%   
%KINSetOptions properties
%(See also the KINSOL User Guide)
% 
%Verbose - verbose output [ {false} | true ]
%   Specifies whether or not KINSOL should output additional information
%MaxNumIter - maximum number of nonlinear iterations [ scalar | {200} ]
%   Specifies the maximum number of iterations that the nonlinar solver is allowed
%   to take.
%FuncRelErr - relative residual error [ scalar | {eps} ]
%   Specifies the realative error in computing f(y) when used in difference
%   quotient approximation of matrix-vector product J(y)*v.
%FuncNormTol - residual stopping criteria [ scalar | {eps^(1/3)} ]
%   Specifies the stopping tolerance on ||fscale*ABS(f(y))||_L-infinity
%ScaledStepTol - step size stopping criteria [ scalar | {eps^(2/3)} ]
%   Specifies the stopping tolerance on the maximum scaled step length:
%                    ||    y_(k+1) - y_k   ||
%                    || ------------------ ||_L-infinity
%                    || |y_(k+1)| + yscale ||
%MaxNewtonStep - maximum Newton step size [ scalar | {0.0} ]
%   Specifies the maximum allowable value of the scaled length of the Newton step.
%InitialSetup - initial call to linear solver setup [ false | {true} ]
%   Specifies whether or not KINSol makes an initial call to the linear solver 
%   setup function.
%MaxNumSetups - [ scalar | {10} ]
%   Specifies the maximum number of nonlinear iterations between calls to the
%   linear solver setup function (i.e. Jacobian/preconditioner evaluation)
%MaxNumSubSetups - [ scalar | {5} ]
%   Specifies the maximum number of nonlinear iterations between checks by the 
%   nonlinear residual monitoring algorithm (specifies length of subintervals).
%   NOTE: MaxNumSetups should be a multiple of MaxNumSubSetups.
%MaxNumBetaFails - maximum number of beta-condition failures [ scalar | {10} ]
%   Specifies the maximum number of beta-condiiton failures in the line search
%   algorithm.
%EtaForm - Inexact Newton method [ Constant | Type2 | {Type1} ]
%   Specifies the method for computing the eta coefficient used in the calculation
%   of the linear solver convergence tolerance (used only if strategy='InexactNEwton'
%   in the call to KINSol):
%      lintol = (eta + eps)*||fscale*f(y)||_L2
%   which is the used to check if the following inequality is satisfied:
%      ||fscale*(f(y)+J(y)*p)||_L2 <= lintol
%   Valid choices are:
%                          | ||f(y_(k+1))||_L2 - ||f(y_k)+J(y_k)*p_k||_L2 |
%   EtaForm='Type1'  eta = ------------------------------------------------
%                                        ||f(y_k)||_L2
%
%                                  [ ||f(y_(k+1))||_L2 ]^alpha
%   EtaForm='Type2'  eta = gamma * [ ----------------- ]
%                                  [  ||f(y_k)||_L2    ]
%   EtaForm='Constant'
%Eta - constant value for eta [ scalar | {0.1} ]
%   Specifies the constant value for eta in the case EtaForm='Constant'.
%EtaAlpha - alpha parameter for eta [ scalar | {2.0} ]
%   Specifies the parameter alpha in the case EtaForm='Type2'
%EtaGamma - gamma parameter for eta [ scalar | {0.9} ]
%   Specifies the parameter gamma in the case EtaForm='Type2'
%MinBoundEps - lower bound on eps [ false | {true} ]
%   Specifies whether or not the value of eps is bounded below by 0.01*FuncNormtol.
%Constraints - solution constraints [ vector ]
%   Specifies additional constraints on the solution components.
%     Constraints(i) =  0 : no constrain on y(i)
%     Constraints(i) =  1 : y(i) >= 0
%     Constraints(i) = -1 : y(i) <= 0
%     Constraints(i) =  2 : y(i) > 0
%     Constraints(i) = -2 : y(i) < 0
%   If Constraints is not specified, no constraints are applied to y.
%
%LinearSolver - Type of linear solver [ {Dense} | Band | GMRES | BiCGStab | TFQMR ]
%   Specifies the type of linear solver to be used for the Newton nonlinear solver. 
%   Valid choices are: Dense (direct, dense Jacobian), GMRES (iterative, scaled 
%   preconditioned GMRES), BiCGStab (iterative, scaled preconditioned stabilized 
%   BiCG), TFQMR (iterative, scaled preconditioned transpose-free QMR).
%   The GMRES, BiCGStab, and TFQMR are matrix-free linear solvers.
%JacobianFn - Jacobian function [ function ]
%   This propeerty is overloaded. Set this value to a function that returns 
%   Jacobian information consistent with the linear solver used (see Linsolver). 
%   If not specified, KINSOL uses difference quotient approximations. 
%   For the Dense linear solver, JacobianFn must be of type KINDenseJacFn and must 
%   return a dense Jacobian matrix. For the iterative linear solvers, GMRES,
%   BiCGStab, or TFQMR, JacobianFn must be of type KINJactimesVecFn and must return 
%   a Jacobian-vector product.
%KrylovMaxDim - Maximum number of Krylov subspace vectors [ scalar | {10} ]
%   Specifies the maximum number of vectors in the Krylov subspace. This property 
%   is used only if an iterative linear solver, GMRES, BiCGStab, or TFQMR is used 
%   (see LinSolver).
%MaxNumRestarts - Maximum number of GMRES restarts [ scalar | {0} ]
%   Specifies the maximum number of times the GMRES (see LinearSolver) solver
%   can be restarted.
%PrecModule - Built-in preconditioner module [ BBDPre | {UserDefined} ]
%   If the PrecModule = 'UserDefined', then the user must provide at least a 
%   preconditioner solve function (see PrecSolveFn)
%   KINSOL provides a built-in preconditioner module, BBDPre which can only be used
%   with parallel vectors. It provide a preconditioner matrix that is block-diagonal 
%   with banded blocks. The blocking corresponds to the distribution of the variable 
%   vector among the processors. Each preconditioner block is generated from the 
%   Jacobian of the local part (on the current processor) of a given function g(t,y) 
%   approximating f(y) (see GlocalFn). The blocks are generated by a difference 
%   quotient scheme on each processor independently. This scheme utilizes an assumed 
%   banded structure with given half-bandwidths, mldq and mudq (specified through 
%   LowerBwidthDQ and UpperBwidthDQ, respectively). However, the banded Jacobian 
%   block kept by the scheme has half-bandwiths ml and mu (specified through 
%   LowerBwidth and UpperBwidth), which may be smaller.
%PrecSetupFn - Preconditioner setup function [ function ]
%   PrecSetupFn specifies an optional function which, together with PrecSolve, 
%   defines a right preconditioner matrix which is an aproximation
%   to the Newton matrix. PrecSetupFn must be of type KINPrecSetupFn. 
%PrecSolveFn - Preconditioner solve function [ function ]
%   PrecSolveFn specifies an optional function which must solve a linear system 
%   Pz = r, for given r. If PrecSolveFn is not defined, the no preconditioning will
%   be used. PrecSolveFn must be of type KINPrecSolveFn.
%GlocalFn - Local right-hand side approximation funciton for BBDPre [ function ]
%   If PrecModule is BBDPre, GlocalFn specifies a required function that
%   evaluates a local approximation to the system function. GlocalFn must
%   be of type KINGlocalFn.
%GcommFn - Inter-process communication function for BBDPre [ function ]
%   If PrecModule is BBDPre, GcommFn specifies an optional function
%   to perform any inter-process communication required for the evaluation of
%   GlocalFn. GcommFn must be of type KINGcommFn.
%LowerBwidth - Jacobian/preconditioner lower bandwidth [ scalar | {0} ]
%   This property is overloaded. If the Band linear solver is used (see LinSolver),
%   it specifies the lower half-bandwidth of the band Jacobian approximation. 
%   If one of the three iterative linear solvers, GMRES, BiCGStab, or TFQMR is used 
%   (see LinSolver) and if the BBDPre preconditioner module in KINSOL is used 
%   (see PrecModule), it specifies the lower half-bandwidth of the retained 
%   banded approximation of the local Jacobian block.
%   LowerBwidth defaults to 0 (no sub-diagonals).
%UpperBwidth - Jacobian/preconditioner upper bandwidth [ scalar | {0} ]
%   This property is overloaded. If the Band linear solver is used (see LinSolver),
%   it specifies the upper half-bandwidth of the band Jacobian approximation. 
%   If one of the three iterative linear solvers, GMRES, BiCGStab, or TFQMR is used 
%   (see LinSolver) and if the BBDPre preconditioner module in KINSOL is used 
%   (see PrecModule), it specifies the upper half-bandwidth of the retained 
%   banded approximation of the local Jacobian block. 
%   UpperBwidth defaults to 0 (no super-diagonals).
%LowerBwidthDQ - BBDPre preconditioner DQ lower bandwidth [ scalar | {0} ]
%   Specifies the lower half-bandwidth used in the difference-quotient Jacobian
%   approximation for the BBDPre preconditioner (see PrecModule).
%UpperBwidthDQ - BBDPre preconditioner DQ upper bandwidth [ scalar | {0} ]
%   Specifies the upper half-bandwidth used in the difference-quotient Jacobian
%   approximation for the BBDPre preconditioner (see PrecModule).
%
%
%   See also
%        KINDenseJacFn, KINJacTimesVecFn
%        KINPrecSetupFn, KINPrecSolveFn
%        KINGlocalFn, KINGcommFn
%

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.3 $Date: 2006/02/02 00:38:58 $

% Based on Matlab's ODESET function

% Print out possible values of properties.
if (nargin == 0) & (nargout == 0)
  fprintf('         Verbose: [ true   | {false} ]\n');
  fprintf('      MaxNumIter: [ scalar | {200} ]\n');
  fprintf('      FuncRelErr: [ scalar | {eps} ]\n');
  fprintf('     FuncNormTol: [ scalar | {eps^(1/3)} ]\n');
  fprintf('   ScaledStepTol: [ scalar | {eps^(2/3)} ]\n');
  fprintf('   MaxNewtonStep: [ scalar | {0.0} ]\n');
  fprintf('    InitialSetup: [ false | {true} ]\n');
  fprintf('    MaxNumSetups: [ scalar | {10} ]\n');
  fprintf(' MaxNumSubSetups: [ scalar | {5} ]\n');  
  fprintf(' MaxNumBetaFails: [ scalar | {10} ]\n');
  fprintf('         EtaForm: [ Constant | Type2 | {Type1} ]\n');
  fprintf('             Eta: [ scalar | {0.1} ]\n');
  fprintf('        EtaAlpha: [ scalar | {2.0} ]\n');
  fprintf('        EtaGamma: [ scalar | {0.9} ]\n');
  fprintf('     MinBoundEps: [ false | {true} ]\n');
  fprintf('     Constraints: [ array of scalar ]\n');
  fprintf('\n');
  fprintf('    LinearSolver: [ {Dense} | Band | GMRES | BiCGStab | TFQMR ]\n');
  fprintf('      JacobianFn: [ function ]\n');
  fprintf('    KrylovMaxDim: [ scalar | {10} ]');
  fprintf('  MaxNumRestarts: [ Classical | {Modified} ]\n');
  fprintf('      PrecModule: [ BBDPre | {UserDefined} ]\n');
  fprintf('     PrecSetupFn: [ function ]\n');
  fprintf('     PrecSolveFn: [ function ]\n');
  fprintf('        GlocalFn: [ function ]\n');
  fprintf('         GcommFn: [ function ]\n');
  fprintf('     LowerBwidth: [ scalar | {0} ]\n');
  fprintf('     UpperBwidth: [ scalar | {0} ]\n');
  fprintf('   LowerBwidthDQ: [ scalar | {0} ]\n');
  fprintf('   UpperBwidthDQ: [ scalar | {0} ]\n');
  fprintf('\n');
  fprintf('\n');
  return;
end

Names = [
    'Verbose         '
    'MaxNumIter      '
    'MaxNumSetups    '
    'MaxNumSubSetups '
    'MaxNumBetaFails '
    'EtaForm         '
    'Eta             '
    'EtaAlpha        '
    'EtaGamma        '
    'MaxNewtonStep   '
    'FuncRelErr      '
    'FuncNormTol     '
    'ScaledStepTol   '
    'InitialSetup    '
    'MinBoundEps     '
    'Constraints     '
    'LinearSolver    '
    'JacobianFn      '
    'PrecType        '
    'PrecModule      '
    'PrecSetupFn     '
    'PrecSolveFn     '
    'GlocalFn        '
    'GcommFn         '
    'KrylovMaxDim    '
    'MaxNumRestarts  '
    'LowerBwidthDQ   '
    'UpperBwidthDQ   '
    'LowerBwidth     '
    'UpperBwidth     '
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
                     'or an options structure\ncreated with KINSetOptions.'], i));
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
