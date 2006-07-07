function [status,t,yB,varargout] = CVodeB(tout,itask)
%CVodeB integrates the backward ODE.
%
%   Usage:  [STATUS, T, YB] = CVodeB ( TOUT, ITASK ) 
%           [STATUS, T, YB, YQB] = CVodeB ( TOUT, ITASK )
%
%   If ITASK is 'Normal', then the solver integrates from its current internal 
%   T value to a point at or beyond TOUT, then interpolates to T = TOUT and returns 
%   YB(TOUT). If ITASK is 'OneStep', then the solver takes one internal time step 
%   and returns in YB the solution at the new internal time. In this case, TOUT 
%   is used only during the first call to CVodeB to determine the direction of 
%   integration and the rough scale of the problem. In either case, the time 
%   reached by the solver is returned in T. 
%
%   If quadratures were computed (see CVodeSetOptions), CVodeB will return their
%   values at T in the vector YQB.
%
%   On return, STATUS is one of the following:
%     0: CVodeB succeeded and no roots were found.
%    -2: One of the inputs to CVodeB is illegal.
%    -4: The solver took mxstep internal steps but could not reach TOUT. 
%        The default value for mxstep is 500.
%   -5:  The solver could not satisfy the accuracy demanded by the user for 
%        some internal step.
%   -6:  Error test failures occurred too many times (MXNEF = 7) during one 
%        internal time step or occurred with |h| = hmin.
%   -7:  Convergence test failures occurred too many times (MXNCF = 10) during 
%        one internal time step or occurred with |h| = hmin.
%   -9:  The linear solver's setup routine failed in an unrecoverable manner.
%  -10:  The linear solver's solve routine failed in an unrecoverable manner.
%  -101: Illegal attempt to call before initializing adjoint sensitivity 
%        (see CVodeMalloc).
%  -104: Illegal attempt to call before CVodeMallocB.
%  -108: Wrong value for TOUT.
%
%   See also CVodeSetOptions, CVodeGetStatsB

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.3 $Date: 2006/03/07 01:19:50 $

mode = 8;
if nargout < 3 | nargout > 4
  disp('CVodeB:: wrong number of arguments');
  return
end

if nargout == 3
  [status,t,yB] = cvm(mode,tout,itask);
elseif nargout == 4
  [status,t,yB,qB] = cvm(mode,tout,itask);
  varargout(1) = {qB};
end