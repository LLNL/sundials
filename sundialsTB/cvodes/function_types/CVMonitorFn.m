%CVMonitorFn - type for user provided monitoring function for forward problems.
%
%   The function MONFUN must be defined as
%       FUNCTION [] = MONFUN(CALL, T, Y, YQ, YS)
%   It is called after every internal CVode step and can be used to
%   monitor the progress of the solver. MONFUN is called with CALL=0
%   from CVodeInit at which time it should initialize itself and it
%   is called with CALL=2 from CVodeFree. Otherwise, CALL=1.
%
%   It receives as arguments the current time T, solution vector Y,
%   and, if they were computed, quadrature vector YQ, and forward 
%   sensitivity matrix YS. If YQ and/or YS were not computed they
%   are empty here.
%
%   If additional data is needed inside MONFUN, it must be defined
%   as
%      FUNCTION NEW_MONDATA = MONFUN(CALL, T, Y, YQ, YS, MONDATA)
%   If the local modifications to the user data structure need to be 
%   saved (e.g. for future calls to MONFUN), then MONFUN must set
%   NEW_MONDATA. Otherwise, it should set NEW_MONDATA=[] 
%   (do not set NEW_MONDATA = DATA as it would lead to unnecessary copying).
%
%   A sample monitoring function, CVodeMonitor, is provided with CVODES.
%
%   See also CVodeSetOptions, CVodeMonitor
%
%   NOTES: 
%   
%   MONFUN is specified through the MonitorFn property in CVodeSetOptions. 
%   If this property is not set, or if it is empty, MONFUN is not used.
%   MONDATA is specified through the MonitorData property in CVodeSetOptions.
%
%   See CVodeMonitor for an implementation example.

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
% $Revision$Date: 2006/03/07 01:19:50 $
