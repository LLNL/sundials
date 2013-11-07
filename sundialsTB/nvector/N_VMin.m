function ret = N_VMin(x,comm)
%N_VMin returns the smallest element of x
%   Usage:  RET = N_VMin ( X [, COMM] )
%
%If COMM is not present, N_VMin returns the minimum value of 
%the local portion of X. Otherwise, it returns the global
%minimum value.

% Radu Serban <radu@llnl.gov>
% LLNS Start Copyright
% Copyright (c) 2013, Lawrence Livermore National Security
% This work was performed under the auspices of the U.S. Department 
% of Energy by Lawrence Livermore National Laboratory in part under 
% Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
% Produced at the Lawrence Livermore National Laboratory.
% All rights reserved.
% For details, see the LICENSE file.
% LLNS End Copyright
% $Revision: 1.1 $Date$

if nargin == 1
  
  ret = min(x);
  
else
  
  lmin = min(x);
  gmin = 0.0;
  MPI_Allreduce(lmin,gmin,'MIN',comm);
  ret = gmin;
  
end