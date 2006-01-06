function [] = KINFree()
%KINFree deallocates memory for the KINSOL solver.
%
%   Usage:  KINFree
%

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$

mode = 6;
kim(mode);
