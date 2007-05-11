% CVODES, an ODE integrator with sensitivity analysis capabilities
%
% The Matlab interface to the SUNDIALS solver CVODES provides access
% to all functionality of the underlying solver, including IVP simulation
% and sensitvity analysis (both forward and adjoint).
%
% The interface consists of several user-callable functions. In addition,
% the user must provide several required and optional user-supplied 
% functions which define the problem to be solved. The user-callable 
% functions and the types of user-supplied functions are listed below.
% For completness, some functions appear more than once.
%
%
% Functions for ODE integration
%
%  CVodeSetOptions     - create an options structure for an ODE problem.
%  CVodeSetQuadOptions - create an options structure for quadrature integration.
%  CVodeInit           - allocate and initialize memory for CVODES.
%  CVodeQuadInit       - allocate and initialize memory for quadrature integration.
%  CVodeReInit         - reinitialize memory for CVODES.
%  CVodeQuadReInit     - reinitialize memory for quadrature integration.
%  CVode               - integrate the ODE problem.
%  CVodeGetStats       - return statistics for the CVODES solver.
%  CVodeGet            - extract data from CVODES memory.
%  CVodeFree           - deallocate memory for the CVODES solver.
%
%
% Functions for forward sensitivity analysis
%
%  CVodeSetOptions     - create an options structure for an ODE problem.
%  CVodeSetQuadOptions - create an options structure for quadrature integration.
%  CVodeSetFSAOptions  - create an options structure for FSA.
%  CVodeInit           - allocate and initialize memory for CVODES.
%  CVodeQuadInit       - allocate and initialize memory for quadrature integration.
%  CVodeSensInit       - allocate and initialize memory for FSA.
%  CVodeReInit         - reinitialize memory for CVODES.
%  CVodeQuadReInit     - reinitialize memory for quadrature integration.
%  CVodeSensReInit     - reinitialize memory for FSA.
%  CVodeSensToggleOff  - temporarily deactivates FSA.
%  CVode               - integrate the ODE problem.
%  CVodeGetStats       - return statistics for the CVODES solver.
%  CVodeGet            - extract data from CVODES memory.
%  CVodeFree           - deallocate memory for the CVODES solver.
%
%
% Functions for adjoint sensitivity analysis
%
%  CVodeSetOptions - create an options structure for an ODE problem.
%  CVodeSetQuadOptions - create an options structure for quadrature integration.
%  CVodeInit           - allocate and initialize memory for the forward problem.
%  CVodeQuadInit       - allocate and initialize memory for forward quadrature integration.
%  CVodeQuadReInit     - reinitialize memory for forward quadrature integration.
%  CVodeReInit         - reinitialize memory for the forward problem.
%  CVodeAdjInit        - allocate and initialize memory for ASA.
%  CVodeInitB          - allocate and initialize a backward problem.
%  CVodeAdjReInit      - reinitialize memory for ASA.
%  CVodeReInitB        - reinitialize a backward problem.
%  CVode               - integrate the forward ODE problem.
%  CVodeB              - integrate the backward problems.
%  CVodeGetStats       - return statistics for the integration of the forward problem.
%  CVodeGetStatsB      - return statistics for the integration of a backward problem.
%  CVodeGet            - extract data from CVODES memory.
%  CVodeFree           - deallocate memory for the CVODES solver.
%
%
% User-supplied function types for forward problems
%
%   CVRhsFn            -  RHS function
%   CVRootFn           -  root-finding function 
%   CVQuadRhsFn        -  quadrature RHS function
%   CVSensRhsFn        -  sensitivity RHS function
%   CVDenseJacFn       -  dense Jacobian function
%   CVBandJacFn        -  banded Jacobian function
%   CVJacTimesVecFn    -  Jacobian times vector function
%   CVPrecSetupFn      -  preconditioner setup function
%   CVPrecSolveFn      -  preconditioner solve function
%   CVGlocalFn         -  RHS approximation function (BBDPre)
%   CVGcomFn           -  communication function (BBDPre)
%   CVMonitorFn        -  monitoring function
%
%
% User-supplied function types for backward problems
%
%   CVRhsFnB           -  RHS function
%   CVQuadRhsFnB       -  quadrature RHS function
%   CVDenseJacFnB      -  dense Jacobian function
%   CVBandJacFnB       -  banded Jacobian function
%   CVJacTimesVecFnB   -  Jacobian times vector function
%   CVPrecSetupFnB     -  preconditioner setup function
%   CVPrecSolveFnB     -  preconditioner solve function
%   CVGlocalFnB        -  RHS approximation function (BBDPre)
%   CVGcomFnB          -  communication function (BBDPre)
%   CVMonitorFnB       -  monitoring function
%
%
% Serial examples provided with the toolbox
%
%   cvdx               - chemical kinetics problem (BDF, Newton, dense)
%   cvbx               - 2D adv-diff problem (BDF, Newton, band)
%   cvkx               - 2D diurnal kinetics problem 
%                        (BDF, Newton, GMRES, user-provided preconditioner)
%   cvkxb              - 2D diurnal kinetics problem 
%                        (BDF, Newton, GMRES, BandPre preconditioner)
%   vdp                - van der Pol problem (BDF, Newton, dense)
%   pleiades           - celestial mechanics problem (Adams, Functional)
%   cvfdx              - FSA for cvdx
%   cvadx              - ASA for cvdx
%   cvhess             - example for multiple backward problems
%
%
% Parallel examples provided with the toolbox
%
%   pvnx               - diagonal ODE example
%   pvfnx              - FSA for 1D adv-diff problem (Adams, Functional)
%   pvkx               - 3D adv-diff with distributed source problem
%                        (BDF, Newton, GMRES, BBDPre preconditioner)
% Use the mpirun function to run any of the parallel examples
%
%
% See also nvector, putils

