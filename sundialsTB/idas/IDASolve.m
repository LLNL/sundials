function [status, t, yy, yp, varargout] = IDASolve(tout,itask)
%IDASolve integrates the DAE.
%
%   Usage: [STATUS, T, YY, YP] = IDASolve ( TOUT, ITASK ) 
%          [STATUS, T, YY, YP, YQ] = IDASolve  (TOUT, ITASK )
%          [STATUS, T, YY, YP, YYS, YPS] = IDASolve ( TOUT, ITASK )
%          [STATUS, T, YY, YP, YQ, YYS, YPS] = IDASolve ( TOUT, ITASK )
%
%   If ITASK is 'Normal', then the solver integrates from its current internal 
%   T value to a point at or beyond TOUT, then interpolates to T = TOUT and 
%   returns YY(TOUT) and YP(TOUT). If ITASK is 'OneStep', then the solver takes 
%   one internal time step and returns in YY and YP the solution at the new 
%   internal time. In this case, TOUT is used only during the first call to 
%   IDASolve to determine the direction of integration and the rough scale of 
%   the problem. In either case, the time reached by the solver is returned in T. 
%   The 'NormalTstop' and 'OneStepTstop' modes are similar to 'Normal' and 
%   'OneStep', respectively, except that the integration never proceeds past 
%   the value tstop.
%
%   If quadratures were computed (see IDASetOptions), IDASolve will return their
%   values at T in the vector YQ.
%
%   If sensitivity calculations were enabled (see IDASetOptions), IDASolve will 
%   return their values at T in the matrix YS.
%
%   On return, STATUS is one of the following:
%     0: IDASolve succeeded and no roots were found.
%     1: IDASolve succeded and returned at tstop.
%     2: IDASolve succeeded, and found one or more roots. 
%    -1: Illegal attempt to call before IDAMalloc
%    -2: One of the inputs to IDASolve is illegal. This includes the situation 
%        when a component of the error weight vectors becomes < 0 during internal 
%        time-stepping.
%    -4: The solver took mxstep internal steps but could not reach TOUT. The 
%        default value for mxstep is 500.
%    -5: The solver could not satisfy the accuracy demanded by the user for some 
%        internal step.
%    -6: Error test failures occurred too many times (MXNEF = 7) during one internal 
%        time step 
%        or occurred with |h| = hmin.
%    -7: Convergence test failures occurred too many times (MXNCF = 10) during one 
%        internal time step or occurred with |h| = hmin.
%    -9: The linear solver's setup routine failed in an unrecoverable manner.
%   -10: The linear solver's solve routine failed in an unrecoverable manner.
%
%
%   See also IDASetOptions, IDAGetStats

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/03/07 01:19:50 $

mode = 10;
if nargout < 4 | nargout > 7
  disp('IDASolve:: wrong number of output arguments');
  return
end

if nargout == 4
  [status,t,yy,yp] = idm(mode,tout,itask);
elseif nargout == 5
  [status,t,yy,yp,yQ] = idm(mode,tout,itask);
  varargout(1) = {yQ};
elseif nargout == 6
  [status,t,yy,yp,yyS,ypS] = idm(mode,tout,itask);
  varargout(1) = {yyS};
  varargout(2) = {ypS};
elseif nargout == 7
  [status,t,yy,yp,yQ,yyS,ypS] = idm(mode,tout,itask);
  varargout(1) = {yQ};
  varargout(2) = {yyS};
  varargout(3) = {ypS};
end