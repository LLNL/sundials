%IDAMonitorFn - type for monitoring function.
%
%   The function MONFUN must be defined as
%      FUNCTION [] = MONFUN(CALL, T, YY, YP, YQ, YYS, YPS)
% 
%   To enable monitoring using a given monitor function MONFUN,
%   use IDASetOptions to set the property 'MonitorFn" to 'MONFUN' 
%   (or to @MONFUN).
%
%   MONFUN is called with the following input arguments:
%
%   o CALL indicates the phase during the integration process at which
%     MONFUN is called:
%     CALL=1 : MONFUN was called at the initial time; this can be either 
%              after IDAInit or after IDAReInit.
%              (typically, MONFUN should perform its own initialization)
%     CALL=2 : MONFUN was called right before a solver reinitializtion.
%              (typically, MONFUN should decide whether to initialize
%              itself or else to continue monitoring)
%     CALL=3 : MONFUN was called during solver finalization.
%              (typically, MONFUN should finalize monitoring)
%     CALL=0 : MONFUN was called after the solver took a successful
%              internal step.
%              (typically, MONFUN should collect and/or display data)
%
%   o T is the current integration time
%
%   o YY and YP are vectors containing the solution and solution 
%     derivative at time T
%
%   o YQ is a vector containing the quadrature variables at time T
%
%   o YYS and YPS are matrices containing the forward sensitivities
%     and their derivatives, respectively, at time T.
%
%   If additional data is needed inside a MONFUN function, then it must 
%   be defined as
%      FUNCTION NEW_MONDATA = MONFUN(CALL, T, YY, YP, YQ, YYS, YPS, MONDATA)
%
%   In this case, the MONFUN function is passed the additional argument
%   MONDATA, the same as that specified through the property 'MonitorData'
%   in IDASetOptions. If the local modifications to the monitor data structure 
%   need to be saved (e.g. for future calls to MONFUN), then MONFUN must set
%   NEW_MONDATA. Otherwise, it should set NEW_MONDATA=[] (do not set 
%   NEW_MONDATA = DATA as it would lead to unnecessary copying).
%
%   NOTES: 
%   
%   1. MONFUN is specified through the MonitorFn property in IDASetOptions. 
%      If this property is not set, or if it is empty, MONFUN is not used.
%      MONDATA is specified through the MonitorData property in IDASetOptions.
%
%   2. If quadrature integration is not enabled, YQ is empty. Similarly, if
%      forward sensitivity analysis is not enabled, YYS and YPS are empty.
%
%   3. When CALL = 2 or 3, all arguments YY, YP, YQ, YYS, and YPS are empty.
%      Moreover, when CALL = 3, T = 0.0
%
%   4. If MONFUN is used on the backward integration phase, YYS and YPS are 
%      always empty.
%
%   See also IDASetOptions, IDAMonitor
%

% Radu Serban <radu@llnl.gov>
% SUNDIALS Copyright Start
% Copyright (c) 2002-2020, Lawrence Livermore National Security
% and Southern Methodist University.
% All rights reserved.
%
% See the top-level LICENSE and NOTICE files for details.
%
% SPDX-License-Identifier: BSD-3-Clause
% SUNDIALS Copyright End
% $Revision$Date: 2007/08/21 17:38:44 $
