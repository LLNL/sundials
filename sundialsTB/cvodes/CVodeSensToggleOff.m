function status = CVodeSensToggleOff()
% CVodeSensToggleOff deactivates sensitivity calculations.
%   It does NOT deallocate sensitivity-related memory so that 
%   sensitivity computations can be later toggled ON (through
%   CVodeSensReInit).
%
%   Usage: CVodeSensToggleOff
%
%   See also: CVodeSensInit, CVodeSensReInit

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
% $Revision$Date: 2007/05/11 18:51:32 $

mode = 18;
status = cvm(mode);
