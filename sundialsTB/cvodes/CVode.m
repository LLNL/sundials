function [status,t,y,varargout] = CVode(tout,itask)
%CVode integrates the ODE.
%
%   Usage: [STATUS, T, Y] = CVode ( TOUT, ITASK ) 
%          [STATUS, T, Y, YS] = CVode ( TOUT, ITASK )
%          [STATUS, T, Y, YQ] = CVode  (TOUT, ITASK )
%          [STATUS, T, Y, YQ, YS] = CVode ( TOUT, ITASK )
%
%   If ITASK is 'Normal', then the solver integrates from its current internal 
%   T value to a point at or beyond TOUT, then interpolates to T = TOUT and returns 
%   Y(TOUT). If ITASK is 'OneStep', then the solver takes one internal time step 
%   and returns in Y the solution at the new internal time. In this case, TOUT 
%   is used only during the first call to CVode to determine the direction of 
%   integration and the rough scale of the problem. In either case, the time 
%   reached by the solver is returned in T. The 'NormalTstop' and 'OneStepTstop' 
%   modes are similar to 'Normal' and 'OneStep', respectively, except that the 
%   integration never proceeds past the value tstop.
%
%   If quadratures were computed (see CVodeSetOptions), CVode will return their
%   values at T in the vector YQ.
%
%   If sensitivity calculations were enabled (see CVodeSetOptions), CVode will 
%   return their values at T in the matrix YS.
%
%   On return, STATUS is one of the following:
%     0: CVode succeeded and no roots were found.
%     1: CVode succeded and returned at tstop.
%     2: CVode succeeded, and found one or more roots. 
%    -1: Illegal attempt to call before CVodeMalloc
%    -2: One of the inputs to CVode is illegal. This includes the situation 
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
%   See also CVodeSetOptions, CVodeGetStats

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.4 $Date: 2006/07/07 19:08:40 $

mode = 20;
if nargout < 3 | nargout > 5
  disp('CVode:: wrong number of output arguments');
  return
end

if nargout == 3
  [status,t,y] = cvm(mode,tout,itask);
elseif nargout == 4
  % v1 can be either yQ or yS
  [status,t,y,v1] = cvm(mode,tout,itask);
  varargout(1) = {v1};
elseif nargout == 5
  [status,t,y,yQ,yS] = cvm(mode,tout,itask);
  varargout(1) = {yQ};
  varargout(2) = {yS};
end