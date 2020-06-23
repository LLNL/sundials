function ret = N_VMaxNorm(x, comm)
%N_VMaxNorm returns the L-infinity norm of x
%
%   Usage:  RET = N_VMaxNorm ( X [, COMM] )
%
%If COMM is not present, N_VMaxNorm returns the L-infinity norm 
%of the local portion of X. Otherwise, it returns the global
%L-infinity norm..

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
  
  ret = norm(x,'inf');
  
else
  
  lnrm = norm(x,'inf');
  gnrm = 0.0;
  MPI_Allreduce(lnrm,gnrm,'MAX',comm);
  ret = gnrm;
  
end