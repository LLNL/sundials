function status = IDASensToggleOff()
% IDASensToggleOff deactivates sensitivity calculations.
%   It does NOT deallocate sensitivity-related memory so that 
%   sensitivity computations can be later toggled ON (through
%   IDASensReInit).
%
%   Usage: IDASensToggleOff
%
%   See also: IDASensInit, IDASensReInit

% Radu Serban <radu@llnl.gov>
% SUNDIALS Copyright Start
% Copyright (c) 2002-2020, Lawrence Livermore National Security
% and Southern Methodist University.
% All rights reserved.
%
% See the top-level LICENSE and NOTICE files for details.
%
% SPDX-License-Identifier: BSD-3-Clause
% SUNDIALS Copyright End
% $Revision$Date: 2007/08/21 17:38:43 $

mode = 18;
status = idm(mode);
