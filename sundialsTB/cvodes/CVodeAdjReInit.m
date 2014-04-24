function status = CVodeAdjReInit()
%CVodeAdjReInit re-initializes memory for ASA with CVODES.
%
%   Usage: CVodeAdjReInit
%

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
% $Revision$Date: 2007/05/11 18:51:31 $

mode = 14;

status = cvm(mode);
