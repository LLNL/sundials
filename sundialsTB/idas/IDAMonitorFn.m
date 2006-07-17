%IDAMonitorFn - type for user provided monitoring function.
%
%   The function MONFUN must be defined as
%       FUNCTION [] = MONFUN(CALL, T, YY, YP, YQ, YYS, YPS)
%   It is called after every internal IDASolve step and can be used to
%   monitor the progress of the solver. MONFUN is called with CALL=0
%   from IDAMalloc at which time it should initialize itself and it
%   is called with CALL=2 from IDAFree. Otherwise, CALL=1.
%
%   It receives as arguments the current time T, solution vectors 
%   YY and YP, and, if they were computed, quadrature vector YQ, and 
%   forward sensitivity matrices YYS and YPS If YQ and/or YYS, YPS were 
%   not computed they are empty here.
%
%   If additional data is needed inside MONFUN, it must be defined
%   as
%      FUNCTION NEW_MONDATA = MONFUN(CALL, T, YY, YP, YQ, YYS, YPS, MONDATA)
%   If the local modifications to the user data structure need to be 
%   saved (e.g. for future calls to MONFUN), then MONFUN must set
%   NEW_MONDATA. Otherwise, it should set NEW_MONDATA=[] 
%   (do not set NEW_MONDATA = DATA as it would lead to unnecessary copying).
%
%   A sample monitoring function, IDAMonitor, is provided with IDAS.
%
%   See also IDASetOptions, IDAMonitor
%
%   NOTES: 
%   
%   MONFUN is specified through the MonitorFn property in IDASetOptions. 
%   If this property is not set, or if it is empty, MONFUN is not used.
%   MONDATA is specified through the MonitorData property in IDASetOptions.
%
%   If MONFUN is used on the backward integration phase, YYS and YPS will 
%   always be empty.
%
%   See IDAMonitor for an example of using MONDATA to write a single
%   monitoring function that works both for the forward and backward
%   integration phases.

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/03/07 01:19:50 $
