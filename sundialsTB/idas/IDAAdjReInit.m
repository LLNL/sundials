function status = IDAAdjReInit()
%IDAAdjReInit re-initializes memory for ASA with CVODES.
%
%   Usage: IDAAdjReInit
%

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
% $Revision: 1.2 $Date: 2007/08/21 17:38:42 $

mode = 14;

status = idm(mode);
