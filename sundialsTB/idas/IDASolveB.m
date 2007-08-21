function [varargout] = IDASolveB(tout,itask)
%IDASolveB integrates the backward DAE.
%
%   Usage:  [STATUS, T, YB] = IDASolveB ( TOUT, ITASK ) 
%           [STATUS, T, YB, YQB] = IDASolveB ( TOUT, ITASK )
%
%   If ITASK is 'Normal', then the solver integrates from its current internal 
%   T value to a point at or beyond TOUT, then interpolates to T = TOUT and returns 
%   YB(TOUT). If ITASK is 'OneStep', then the solver takes one internal time step 
%   and returns in YB the solution at the new internal time. In this case, TOUT 
%   is used only during the first call to CVodeB to determine the direction of 
%   integration and the rough scale of the problem. In either case, the time 
%   reached by the solver is returned in T. 
%
%   If quadratures were computed (see CVodeQuadInitB), CVodeB will return their
%   values at T in the vector YQB.
%
%   In ITASK =' Normal' mode, to obtain solutions at specific times T0,T1,...,TFINAL
%   (all increasing or all decreasing) use TOUT = [T0 T1  ... TFINAL]. In this case
%   the output arguments YB and YQB are matrices, each column representing the solution
%   vector at the corresponding time returned in the vector T.
%
%   If more than one backward problem was defined, the return arguments are cell
%   arrays, with T{IDXB}, YB{IDXB}, and YQB{IDXB} corresponding to the backward
%   problem with index IDXB (as returned by CVodeInitB).
%
%   On return, STATUS is one of the following:
%     0: IDASolveB succeeded.
%     1: IDASolveB succeded and return at a tstop value (internally set).
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
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.3 $Date: 2007/02/05 20:23:47 $

mode = 21;

if nargin ~= 2
  error('Wrong number of input arguments');
end

if nargout < 3 || nargout > 4
  error('Wrong number of output arguments');
end

varargout = cell (nargout, 1);
[varargout{:}] = idm(mode,tout,itask);

