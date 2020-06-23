function ret = N_VMax(x,comm)
%N_VMax returns the largest element of x
%
%   Usage:  RET = N_VMax ( X [, COMM] )
%
%If COMM is not present, N_VMax returns the maximum value of 
%the local portion of X. Otherwise, it returns the global
%maximum value.

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
% $Revision$Date$

if nargin == 1
  
  ret = max(x);
  
else
  
  lmax = max(x);
  gmax = 0.0;
  MPI_Allreduce(lmax,gmax,'MAX',comm);
  ret = gmax;
  
end