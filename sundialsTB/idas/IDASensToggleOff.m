function [] = IDASensToggleOff()
% IDASensToggleOff deactivates sensitivity calculations.
%   It does NOT deallocate sensitivity-related memory so that 
%   sensitivity computations can be later toggled ON (through
%   IDASensReInit).
%
%   Usage: IDASensToggleOff
%
%   See also: IDASensInit, IDASensReInit

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2007/02/05 20:23:47 $

mode = 18;
idm(mode);
