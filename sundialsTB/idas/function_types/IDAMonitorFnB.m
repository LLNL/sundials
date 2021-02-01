%IDAMonitorFnB - type of monitoring function for backward problems.
%
%   The function MONFUNB must be defined as
%       FUNCTION [] = MONFUNB(CALL, IDXB, T, Y, YQ)
%   It is called after every internal IDASolveB step and can be used to
%   monitor the progress of the solver. MONFUNB is called with CALL=0
%   from IDAInitB at which time it should initialize itself and it
%   is called with CALL=2 from IDAFree. Otherwise, CALL=1.
%
%   It receives as arguments the index of the backward problem (as
%   returned by IDAInitB), the current time T, solution vector Y,
%   and, if it was computed, the quadrature vector YQ. If quadratures
%   were not computed for this backward problem, YQ is empty here.
%
%   If additional data is needed inside MONFUNB, it must be defined
%   as
%      FUNCTION NEW_MONDATA = MONFUNB(CALL, IDXB, T, Y, YQ, MONDATA)
%   If the local modifications to the user data structure need to be 
%   saved (e.g. for future calls to MONFUNB), then MONFUNB must set
%   NEW_MONDATA. Otherwise, it should set NEW_MONDATA=[] 
%   (do not set NEW_MONDATA = DATA as it would lead to unnecessary copying).
%
%   A sample monitoring function, IDAMonitorB, is provided with CVODES.
%
%   See also IDASetOptions, IDAMonitorB
%
%   NOTES: 
%   
%   MONFUNB is specified through the MonitorFn property in IDASetOptions. 
%   If this property is not set, or if it is empty, MONFUNB is not used.
%   MONDATA is specified through the MonitorData property in IDASetOptions.
%
%   See IDAMonitorB for an implementation example.

% Radu Serban <radu@llnl.gov>
% SUNDIALS Copyright Start
% Copyright (c) 2002-2021, Lawrence Livermore National Security
% and Southern Methodist University.
% All rights reserved.
%
% See the top-level LICENSE and NOTICE files for details.
%
% SPDX-License-Identifier: BSD-3-Clause
% SUNDIALS Copyright End
% $Revision$Date: 2007/05/11 18:51:33 $
