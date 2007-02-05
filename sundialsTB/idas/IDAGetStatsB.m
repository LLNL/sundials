function si = IDAGetStatsB()
%IDAGetStatsB returns run statistics for the backward IDAS solver.
%
%   Usage: STATS = IDAGetStatsB
%
%Fields in the structure STATS
%
%o nst - number of integration steps
%o nre - number of residual function evaluations
%o nsetups - number of linear solver setup calls
%o netf - number of error test failures
%o nni - number of nonlinear solver iterations
%o ncfn - number of convergence test failures
%o qlast - last method order used
%o qcur - current method order
%o h0used - actual initial step size used
%o hlast - last step size used
%o hcur - current step size
%o tcur - current time reached by the integrator
%o QuadInfo - structure with quadrature integration statistics
%o LSInfo - structure with linear solver statistics
%
%The structure LSinfo has different fields, depending on the linear solver used.
%
%If quadratures were present, the structure QuadInfo has the following fields
%
%o nfQe - number of quadrature integrand function evaluations
%o netfQ - number of error test failures for quadrature variables
%
%  Fields in LSinfo for the 'Dense' linear solver
%
%o name - 'Dense'
%o njeD - number of Jacobian evaluations
%o nreD - number of residual function evaluations for difference-quotient
%         Jacobian approximation
%
%  Fields in LSinfo for the 'Band' linear solver
%
%o name - 'Band'
%o njeB - number of Jacobian evaluations
%o nreB - number of residual function evaluations for difference-quotient
%         Jacobian approximation
%
%  Fields in LSinfo for the 'GMRES' and 'BiCGStab' linear solvers
%
%o name - 'GMRES' or 'BiCGStab'
%o nli - number of linear solver iterations
%o npe - number of preconditioner setups
%o nps - number of preconditioner solve function calls
%o ncfl - number of linear system convergence test failures
%o njeSG - number of Jacobian-vector product evaluations
%o nreSG -  number of residual function evaluations for difference-quotient
%          Jacobian-vector product approximation

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/07/17 16:49:50 $

mode = 31;
si = idm(mode);
