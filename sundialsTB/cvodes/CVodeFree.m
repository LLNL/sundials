function CVodeFree()
%CVodeFree deallocates memory for the CVODES solver.
%
%   Usage:  CVodeFree

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
% $Revision$Date: 2006/10/11 18:12:36 $

mode = 40;
cvm(mode);
