function status = CVodeAdjReInit()
%CVodeAdjReInit re-initializes memory for ASA with CVODES.
%
%   Usage: CVodeAdjReInit
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
% $Revision$Date: 2007/05/11 18:51:31 $

mode = 14;

status = cvm(mode);
