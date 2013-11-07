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

