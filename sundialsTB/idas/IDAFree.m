function [] = IDAFree()
%IDAFree deallocates memory for the IDAS solver.
%
%   Usage:  IDAFree

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/03/07 01:19:50 $

mode = 30;
idm(mode);
