function [] = IDAFree()
%IDAFree deallocates memory for the IDAS solver.
%
%   Usage:  IDAFree

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.3 $Date: 2007/02/05 20:23:46 $

mode = 40;
idm(mode);
