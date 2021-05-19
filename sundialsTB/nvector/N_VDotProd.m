function ret = N_VDotProd(x,y,comm)
%N_VDotProd returns the dot product of two vectors
%
%   Usage:  RET = N_VDotProd ( X, Y [, COMM] )
%
%If COMM is not present, N_VDotProd returns the dot product of the
%local portions of X and Y. Otherwise, it returns the global dot
%product.

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
  
  ret = dot(x,y);
  
else
  
  ldot = dot(x,y);
  gdot = 0.0;
  MPI_Allreduce(ldot,gdot,'SUM',comm);
  ret = gdot;
  
end