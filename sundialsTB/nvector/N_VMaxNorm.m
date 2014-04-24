function ret = N_VMaxNorm(x, comm)
%N_VMaxNorm returns the L-infinity norm of x
%
%   Usage:  RET = N_VMaxNorm ( X [, COMM] )
%
%If COMM is not present, N_VMaxNorm returns the L-infinity norm 
%of the local portion of X. Otherwise, it returns the global
%L-infinity norm..

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
  
  ret = norm(x,'inf');
  
else
  
  lnrm = norm(x,'inf');
  gnrm = 0.0;
  MPI_Allreduce(lnrm,gnrm,'MAX',comm);
  ret = gnrm;
  
end