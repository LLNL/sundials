function CVodeFree()
%CVodeFree deallocates memory for the CVODES solver.
%
%   Usage:  CVodeFree

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.4 $Date: 2006/10/11 18:12:36 $

mode = 40;
cvm(mode);
