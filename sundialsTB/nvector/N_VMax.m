function ret = N_VMax(x,comm)
%N_VMax returns the largest element of x
%
%   Usage:  RET = N_VMax ( X [, COMM] )
%
%If COMM is not present, N_VMax returns the maximum value of 
%the local portion of X. Otherwise, it returns the global
%maximum value.

% Radu Serban <radu@llnl.gov>
% LLNS Copyright Start
% Copyright (c) 2014, Lawrence Livermore National Security
% This work was performed under the auspices of the U.S. Department 
% of Energy by Lawrence Livermore National Laboratory in part under 
% Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
% Produced at the Lawrence Livermore National Laboratory.
% All rights reserved.
% For details, see the LICENSE file.
% LLNS Copyright End
% $Revision$Date$

if nargin == 1
  
  ret = max(x);
  
else
  
  lmax = max(x);
  gmax = 0.0;
  MPI_Allreduce(lmax,gmax,'MAX',comm);
  ret = gmax;
  
end