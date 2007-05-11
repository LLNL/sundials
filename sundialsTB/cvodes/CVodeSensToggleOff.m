function [] = CVodeSensToggleOff()
% CVodeSensToggleOff deactivates sensitivity calculations.
%   It does NOT deallocate sensitivity-related memory so that 
%   sensitivity computations can be later toggled ON (through
%   CVodeSensReInit).
%
%   Usage: CVodeSensToggleOff
%
%   See also: CVodeSensInit, CVodeSensReInit

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2007, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/11/25 19:57:25 $

mode = 18;
cvm(mode);
