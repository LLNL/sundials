% IDAS, a DAE integrator with sensitivity analysis capabilities
%
% NOTE:
%    FUNCTIONS MARKED WITH * ARE NOT CURRENTLY AVAILABLE
%
% Functions
%   IDASetOptions - creates an options structure for IDAS.
%   IDAMalloc     - allocates and initializes memory for IDAS.
% *  IDAMallocB    - allocates and initializes backward memory for IDAS.
%   IDASolve      - integrates the DAE
% *  IDASolveB     - integrates the backward DAE.
%   IDAGetStats   - returns statistics for the IDAS solver
% *  IDAGetStatsB  - returns statistics for the backward IDAS solver
% *  IDAGet        - extracts data from IDAS memory
%   IDAFree       - deallocates memory for the IDAS solver.
%
% User-supplied function types   
%   IDAResFn         -  residual function
%   IDARootFn        -  root-finding function 
% *  IDAQuadRhsFn     -  quadrature RHS function
%   IDADenseJacFn    -  dense Jacobian function
%   IDABandJacFn     -  banded Jacobian function
%   IDAJacTimesVecFn -  Jacobian times vector function
%   IDAPrecSetupFn   -  preconditioner setup function
%   IDAPrecSolveFn   -  preconditioner solve function
%   IDAGlocalFn      -  residual approximation function (BBDPre)
%   IDAGcomFn        -  communication function (BBDPre)
% *  IDASensResFn     -  sensitivity residual function
%   IDAMonitorFn     -  monitoring function
%
% Serial examples
%   idadenx  - Robertson chemical kinetics example
%   idabanx  - 1-D heat equation
%
% Parallel examples
%   N/A
%
% Use the mpirun function to run any of the parallel examples
%
% See also nvector, putils

