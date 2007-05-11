function CVodeAdjReInit()
%CVodeAdjReInit re-initializes memory for ASA with CVODES.
%
%   Usage: CVodeAdjReInit
%

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/03/15 19:31:25 $

mode = 14;

cvm(mode);
