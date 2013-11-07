function ret = N_VL1Norm(x,comm)
%N_VL1Norm returns the L1 norm of x
%
%   Usage:  RET = N_VL1Norm ( X [, COMM] )
%
%If COMM is not present, N_VL1Norm returns the L1 norm of 
%the local portion of X. Otherwise, it returns the global
%L1 norm..

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
  
  ret = norm(x,1);
  
else
  
  lnrm = norm(x,1);
  gnrm = 0.0;
  MPI_Allreduce(lnrm,gnrm,'MAX',comm);
  ret = gnrm;
  
end