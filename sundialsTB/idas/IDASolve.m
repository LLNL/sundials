function [varargout] = IDASolve(tout,itask)
%IDASolve integrates the DAE.
%
%   Usage: [STATUS, T, Y] = IDASolve ( TOUT, ITASK ) 
%          [STATUS, T, Y, YQ] = IDASolve  (TOUT, ITASK )
%          [STATUS, T, Y, YS] = IDASolve ( TOUT, ITASK )
%          [STATUS, T, Y, YQ, YS] = IDASolve ( TOUT, ITASK )
%
%   If ITASK is 'Normal', then the solver integrates from its current internal 
%   T value to a point at or beyond TOUT, then interpolates to T = TOUT and returns 
%   Y(TOUT). If ITASK is 'OneStep', then the solver takes one internal time step 
%   and returns in Y the solution at the new internal time. In this case, TOUT 
%   is used only during the first call to IDASolve to determine the direction of 
%   integration and the rough scale of the problem. In either case, the time 
%   reached by the solver is returned in T.
%
%   If quadratures were computed (see IDAQuadInit), IDASolve will return their
%   values at T in the vector YQ.
%
%   If sensitivity calculations were enabled (see IDASensInit), IDASolve will 
%   return their values at T in the matrix YS. Each row in the matrix YS
%   represents the sensitivity vector with respect to one of the problem parameters.
%
%   In ITASK =' Normal' mode, to obtain solutions at specific times T0,T1,...,TFINAL
%   (all increasing or all decreasing) use TOUT = [T0 T1  ... TFINAL]. In this case
%   the output arguments Y and YQ are matrices, each column representing the solution
%   vector at the corresponding time returned in the vector T. If computed, the 
%   sensitivities are eturned in the 3-dimensional array YS, with YS(:,:,I) representing
%   the sensitivity vectors at the time T(I).
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
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.3 $Date: 2007/02/05 20:23:47 $


mode = 20;

if nargin ~= 2
  error('Wrong number of input arguments');
end

if nargout < 3 || nargout > 5
  error('Wrong number of output arguments');
end

varargout = cell (nargout, 1);

[varargout{:}] = idm(mode,tout,itask);
