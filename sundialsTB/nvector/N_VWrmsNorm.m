function ret = N_VWrmsNorm(x,w,comm)
%N_VWrmsNorm returns the weighted root mean square norm of x 
%with weight vector w: 
%   sqrt [(sum (i = 0 to N-1) {(x[i]*w[i])^2})/N]
%
%   Usage:  RET = N_VWrmsNorm ( X, W [, COMM] )   
%
%If COMM is not present, N_VWrmsNorm returns the WRMS norm
%of the local portion of X. Otherwise, it returns the global
%WRMS norm..

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

if nargin == 2
  
  ret = dot(x.^2,w.^2);
  ret = sqrt(ret/length(x));
  
else
  
  lnrm = dot(x.^2,w.^2);
  gnrm = 0.0;
  MPI_Allreduce(lnrm,gnrm,'SUM',comm);

  ln = length(x);
  gn = 0;
  MPI_Allreduce(ln,gn,'SUM',comm);
  
  ret = sqrt(gnrm/gn);
  
end

