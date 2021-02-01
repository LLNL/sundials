function ret = N_VWL2Norm(x,w,comm)
%N_VWL2Norm returns the weighted Euclidean L2 norm of x 
%   with weight vector w:
%   sqrt [(sum (i = 0 to N-1) {(x[i]*w[i])^2})]
%
%   Usage:  RET = N_VWL2Norm ( X, W [, COMM] )
%
%If COMM is not present, N_VWL2Norm returns the weighted L2
%norm of the local portion of X. Otherwise, it returns the 
%global weighted L2 norm..

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

if nargin == 2
  
  ret = dot(x.^2,w.^2);
  ret = sqrt(ret);
  
else
  
  lnrm = dot(x.^2,w.^2);
  gnrm = 0.0;
  MPI_Allreduce(lnrm,gnrm,'SUM',comm);
  
  ret = sqrt(gnrm);
  
end