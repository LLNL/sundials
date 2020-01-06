%IDAGcommFnB - type for communication function (BBDPre) for backward problems.
%
%   The function GCOMFUNB must be defined either as
%        FUNCTION FLAG = GCOMFUNB(T, YY, YP, YYB, YPB)
%   or as
%        FUNCTION [FLAG, NEW_DATA] = GCOMFUNB(T, YY, YP, YYB, YPB, DATA)
%   depending on whether a user data structure DATA was specified in
%   IDAInit. 
%
%   The function GCOMFUNB must set FLAG=0 if successful, FLAG<0 if an
%   unrecoverable failure occurred, or FLAG>0 if a recoverable error
%   occurred.
%
%   See also IDAGlocalFnB, IDAGcommFn, IDASetOptions
%
%   NOTES:
%     GCOMFUNB is specified through the GcommFn property in IDASetOptions 
%     and is used only if the property PrecModule is set to 'BBDPre'.
%
%     Each call to GCOMFUNB is preceded by a call to the residual function
%     DAEFUN with the same arguments T, YY, YP and YYB and YPB.
%     Thus GCOMFUNB can omit any communication done by DAEFUNB if relevant
%     to the evaluation of G by GLOCFUNB. If all necessary communication 
%     was done by DAEFUNB, GCOMFUNB need not be provided.     

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
% $Revision$Date: 2007/08/21 17:38:44 $
