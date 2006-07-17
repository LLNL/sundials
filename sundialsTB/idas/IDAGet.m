function varargout = IDAGet(key, varargin)
%IDAGet extracts data from the IDAS solver memory.
%
%   Usage: RET = IDAGet ( KEY [, P1 [, P2] ... ]) 
%
%   IDAGet returns internal IDAS information based on KEY. For some values
%   of KEY, additional arguments may be required and/or more than one output is
%   returned.
%
%   KEY is a string and should be one of:
%    o ErrorWeights - Returns a vector containing the current error weights.
%       EWT = IDAGet('ErrorWeights')
%    o CheckPointsInfo - Returns an array of structures with check point information.
%       CK = IDAGet('CheckPointInfo)
%    o CurrentCheckPoint - Returns the address of the active check point
%       ADDR = IDAGet('CurrentCheckPoint');
%    o DataPointInfo - Returns information stored for interpolation at the I-th data
%       point in between the current check points. The index I must be passed through
%       the agument P1.
%       If the interpolation type was Hermite (see IDASetOptions), it returns two
%       vectors, Y and YD:
%       [Y, YD] = IDAGet('DataPointInfo', I)

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date: 2006/03/07 01:19:50 $

mode = 22;

if strcmp(key, 'ErrorWeights')
  ewt = idm(mode,2);
  varargout(1) = {ewt};
elseif strcmp(key, 'CheckPointsInfo')
  ck = idm(mode,4);
  varargout(1) = {ck};
elseif strcmp(key, 'CurrentCheckPoint')
  addr = idm(mode, 5);
  varargout(1) = {addr};
elseif strcmp(key, 'DataPointInfo')
  i = varargin{1};
  varargout = idm(mode,6,i);
else
  error('IDAGet:: Unrecognized key');
end