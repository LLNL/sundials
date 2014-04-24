function status = IDAReInit(t0,yy0,yp0,options)
%IDAReInit reinitializes memory for IDAS.
%   where a prior call to IDAInit has been made with the same
%   problem size N. IDAReInit performs the same input checking
%   and initializations that IDAInit does, but it does no 
%   memory allocation, assuming that the existing internal memory 
%   is sufficient for the new problem.
%
%   Usage: IDAReInit ( T0, YY0, YP0 [, OPTIONS ] ) 
%
%   T0       is the initial value of t.
%   YY0      is the initial condition vector y(t0).  
%   YP0      is the initial condition vector y'(t0).  
%   OPTIONS  is an (optional) set of integration options, created with
%            the IDASetOptions function. 
%
%  See also: IDASetOptions, IDAInit

% Radu Serban <radu@llnl.gov>
% LLNS Copyright Start
% Copyright (c) 2014, Lawrence Livermore National Security
% This work was performed under the auspices of the U.S. Department 
% of Energy by Lawrence Livermore National Laboratory in part under 
% Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
% Produced at the Lawrence Livermore National Laboratory.
% All rights reserved.
% For details, see the LICENSE file.
% LLNS Copyright End
% $Revision$Date: 2007/08/21 17:38:42 $

mode = 11;

if nargin < 3
  error('Too few input arguments');
end

if nargin < 4
  options = [];
end

status = idm(mode, t0, yy0, yp0, options);
