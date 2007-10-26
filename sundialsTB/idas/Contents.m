% IDAS, a DAE integrator with sensitivity analysis capabilities
%
% The Matlab interface to the SUNDIALS solver IDAS provides access
% to all functionality of the underlying solver, including IVP simulation
% and sensitvity analysis (both forward and adjoint).
%
% The interface consists of several user-callable functions. In addition,
% the user must provide several required and optional user-supplied 
% functions which define the problem to be solved. The user-callable 
% functions and the types of user-supplied functions are listed below.
% For completness, some functions appear more than once.
%
% Functions for DAE integration
%
%  IDASetOptions     - create an options structure for IDAS.
%  IDASetQuadOptions - create an options structure for quadrature integration.
%  IDAInit           - allocate and initialize memory for IDAS.
%  IDAQuadInit       - allocate and initialize memory for quadrature integration.
%  IDAReInit         - reinitialize memory for IDAS.
%  IDAQuadReInit     - reinitialize memory for quadrature integration.
%  IDACalcIC         - compute consistent initial conditions.
%  IDASolve          - integrate the DAE problem.
%  IDAGetStats       - return statistics for the IDAS solver.
%  IDAGet            - extract data from IDAS memory.
%  IDAFree           - deallocate memory for the IDAS solver.
%
% Functions for forward sensitivity analysis
%
%  IDASetOptions     - create an options structure for an DAE problem.
%  IDAQuadSetOptions - create an options structure for quadrature integration.
%  IDASensSetOptions - create an options structure for FSA.
%  IDAInit           - allocate and initialize memory for IDAS.
%  IDAQuadInit       - allocate and initialize memory for quadrature integration.
%  IDASensInit       - allocate and initialize memory for FSA.
%  IDAReInit         - reinitialize memory for IDAS.
%  IDAQuadReInit     - reinitialize memory for quadrature integration.
%  IDASensReInit     - reinitialize memory for FSA.
%  IDASensToggleOff  - temporarily deactivates FSA.
%  IDASetIC          - compute consistent initial conditions.
%  IDASolve          - integrate the DAE problem.
%  IDAGetStats       - return statistics for the IDAS solver.
%  IDAGet            - extract data from IDAS memory.
%  IDAFree           - deallocate memory for the IDAS solver.
%
% Functions for adjoint sensitivity analysis
%
%  IDASetOptions     - create an options structure for an DAE problem.
%  IDAQuadSetOptions - create an options structure for quadrature integration.
%  IDAInit           - allocate and initialize memory for the forward problem.
%  IDAQuadInit       - allocate and initialize memory for forward quadrature integration.
%  IDAQuadReInit     - reinitialize memory for forward quadrature integration.
%  IDAReInit         - reinitialize memory for the forward problem.
%  IDAAdjInit        - allocate and initialize memory for ASA.
%  IDAInitB          - allocate and initialize a backward problem.
%  IDAAdjReInit      - reinitialize memory for ASA.
%  IDAReInitB        - reinitialize a backward problem.
%  IDASetIC          - compute consistent initial conditions.
%  IDASetICb         - compute consistent final conditions for backward problem.
%  IDASolve          - integrate the forward DAE problem.
%  IDASolveB         - integrate the backward problems.
%  IDAGetStats       - return statistics for the integration of the forward problem.
%  IDAGetStatsB      - return statistics for the integration of a backward problem.
%  IDAGet            - extract data from IDAS memory.
%  IDAFree           - deallocate memory for the IDAS solver.
%
% User-supplied function types for forward problems
%
%   IDAResFn            -  DAE residual function
%   IDARootFn           -  root-finding function 
%   IDAQuadRhsFn        -  quadrature RHS function
%   IDASensResFn        -  sensitivity DAE residual function
%   IDADenseJacFn       -  dense Jacobian function
%   IDABandJacFn        -  banded Jacobian function
%   IDAJacTimesVecFn    -  Jacobian times vector function
%   IDAPrecSetupFn      -  preconditioner setup function
%   IDAPrecSolveFn      -  preconditioner solve function
%   IDAGlocalFn         -  RHS approximation function (BBDPre)
%   IDAGcomFn           -  communication function (BBDPre)
%   IDAMonitorFn        -  monitoring function
%
% User-supplied function types for backward problems
%
%   IDAResFnB           -  backard DAE residual function
%   IDAQuadRhsFnB       -  quadrature RHS function
%   IDADenseJacFnB      -  dense Jacobian function
%   IDABandJacFnB       -  banded Jacobian function
%   IDAJacTimesVecFnB   -  Jacobian times vector function
%   IDAPrecSetupFnB     -  preconditioner setup function
%   IDAPrecSolveFnB     -  preconditioner solve function
%   IDAGlocalFnB        -  RHS approximation function (BBDPre)
%   IDAGcomFnB          -  communication function (BBDPre)
%   IDAMonitorFnB       -  monitoring function
%
% Serial examples provided with the toolbox
%
%   midasRoberts_dns      -  chemical kinetics problem (index-1 DAE)
%   midasRoberts_ASAi_dns -  ASA for the robertson problem
%   midasBruss_dns        -  2D, 2-species, time dependent PDE (index-1 DAE)
%   midasBruss_ASA_dns    -  ASA for the brusselator example
%   midasHeat2D_bnd       -  2D heat problem
%   midasPendI1_dns       -  simple pendulum example (index-1 DAE)
%   midasPendI2_dns       -  simple pendulum example (stabilized index-2 DAE)
%   midasSlCrank_dns      -  slider-crank example (stabilized index-2 DAE)
%   midasSlCrank_FSA_dns  -  FSA for the slider-crank example
%   midasReInit_dns       -  integration over solution discontinuities
%
% Parallel examples provided with the toolbox
%
%    N/A
% Use the mpirun function to run any of the parallel examples
%
% See also nvector, putils

