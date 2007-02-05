function [status,t,yyB,ypB,varargout] = IDASolveB(tout,itask)
%IDASolveB integrates the backward DAE.
%
%   Usage:  [STATUS, T, YYB, YPB] = IDASolveB ( TOUT, ITASK ) 
%           [STATUS, T, YYB, YPB, YQB] = IDASolveB ( TOUT, ITASK )
%
%   If ITASK is 'Normal', then the solver integrates from its current internal 
%   T value to a point at or beyond TOUT, then interpolates to T = TOUT and 
%   returns YYB(TOUT) and YPB(TOUT). If ITASK is 'OneStep', then the solver 
%   takes one internal time step and returns in YYB and YPB the solution at 
%   the new internal time. In this case, TOUT is used only during the first 
%   call to IDASolveB to determine the direction of integration and the rough 
%   scale of the problem. In either case, the time reached by the solver is 
%   returned in T. 
%
%   If quadratures were computed (see IDASetOptions), IDASolveB will return their
%   values at T in the vector YQB.
%
%   On return, STATUS is one of the following:
%     0: IDASolveB succeeded and no roots were found.
%    -2: One of the inputs to IDASolveB is illegal.
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
%        (see IDAMalloc).
%  -104: Illegal attempt to call before IDAMallocB.
%  -108: Wrong value for TOUT.
%
%   See also IDASetOptions, IDAGetStatsB

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/07/17 16:49:50 $

mode = 23;
if nargout < 4 | nargout > 5
  disp('IDASolveB:: wrong number of arguments');
  return
end

if nargout == 4
  [status,t,yyB,ypB] = idm(mode,tout,itask);
elseif nargout == 4
  [status,t,yyB,ypB,yqB] = idm(mode,tout,itask);
  varargout(1) = {yqB};
end