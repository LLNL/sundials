function KINFree()
%KINFree deallocates memory for the KINSOL solver.
%
%   Usage:  KINFree

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
% $Revision$Date: 2006/01/06 19:00:02 $

mode = 6;
kim(mode);
