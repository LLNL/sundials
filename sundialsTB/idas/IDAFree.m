function [] = IDAFree()
%IDAFree deallocates memory for the IDAS solver.
%
%   Usage:  IDAFree

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
% $Revision$Date: 2007/02/05 20:23:46 $

mode = 40;
idm(mode);
