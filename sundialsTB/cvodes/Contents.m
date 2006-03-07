% CVODES, an ODE integrator with sensitivity analysis capabilities
%
% Functions
%   CVodeSetOptions - creates an options structure for CVODES.
%   CVodeMalloc     - allocates and initializes memory for CVODES.
%   CVodeMallocB    - allocates and initializes backward memory for CVODES.
%   CVode           - integrates the ODE
%   CVodeB          - integrates the backward ODE.
%   CVodeGetStats   - returns statistics for the CVODES solver
%   CVodeGetStatsB  - returns statistics for the backward CVODES solver
%   CVodeGet        - extracts data from CVODES memory
%   CVodeFree       - deallocates memory for the CVODES solver.
%
% User-supplied function types   
%   CVRhsFn         -  RHS function
%   CVRootFn        -  root-finding function 
%   CVQuadRhsFn     -  quadrature RHS function
%   CVDenseJacFn    -  dense Jacobian function
%   CVBandJacFn     -  banded Jacobian function
%   CVJacTimesVecFn -  Jacobian times vector function
%   CVPrecSetupFn   -  preconditioner setup function
%   CVPrecSolveFn   -  preconditioner solve function
%   CVGlocalFn      -  RHS approximation function (BBDPre)
%   CVGcomFn        -  communication function (BBDPre)
%   CVSensRhsFn     -  sensitivity RHS function
%   CVMonitorFn     -  monitoring function
%
% Serial examples
%   cvdx     - chemical kinetics problem (BDF, Newton, dense)
%   cvbx     - 2D adv-diff problem (BDF, Newton, band)
%   cvkx     - 2D diurnal kinetics problem 
%              (BDF, Newton, GMRES, user-provided preconditioner)
%   cvkxb    - 2D diurnal kinetics problem 
%              (BDF, Newton, GMRES, BandPre preconditioner)
%   vdp      - van der Pol problem (BDF, Newton, dense)
%   pleiades - celestial mechanics problem (Adams, Functional)
%   cvfdx    - FSA for cvdx
%   cvadx    - ASA for cvdx
%
% Parallel examples
%   pvnx  - diagonal ODE example
%   pvfnx - FSA for 1D adv-diff problem (Adams, Functional)
%   pvkx  - 3D adv-diff with distributed source problem
%           (BDF, Newton, GMRES, BBDPre preconditioner)
% Use the mpirun function to run any of the parallel examples
%
% See also nvector, putils

