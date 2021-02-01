function ret = N_VL1Norm(x,comm)
%N_VL1Norm returns the L1 norm of x
%
%   Usage:  RET = N_VL1Norm ( X [, COMM] )
%
%If COMM is not present, N_VL1Norm returns the L1 norm of 
%the local portion of X. Otherwise, it returns the global
%L1 norm..

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
% $Revision$Date$

if nargin == 1
  
  ret = norm(x,1);
  
else
  
  lnrm = norm(x,1);
  gnrm = 0.0;
  MPI_Allreduce(lnrm,gnrm,'MAX',comm);
  ret = gnrm;
  
end