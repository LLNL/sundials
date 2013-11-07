function status = KINInit(fct, n, options)
%KINInit allocates and initializes memory for KINSOL.
%
%   Usage:   KINInit ( SYSFUN, N [, OPTIONS ] );
%
%   SYSFUN   is a function defining the nonlinear problem f(y) = 0.
%            This function must return a column vector FY containing the
%            current value of the residual
%   N        is the (local) problem dimension.
%   OPTIONS  is an (optional) set of integration options, created with
%            the KINSetOptions function. 
%
%   See also: KINSetOptions, KINSysFn 

% Radu Serban <radu@llnl.gov>
% LLNS Start Copyright
% Copyright (c) 2013, Lawrence Livermore National Security
% This work was performed under the auspices of the U.S. Department 
% of Energy by Lawrence Livermore National Laboratory in part under 
% Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
% Produced at the Lawrence Livermore National Laboratory.
% All rights reserved.
% For details, see the LICENSE file.
% LLNS End Copyright
% $Revision: 1.1 $Date: 2006/01/06 19:00:02 $

mode = 1;

if nargin < 2
  error('Too few input arguments');
end

if nargin < 3
  options = [];
end

status = kim(mode, fct, n, options);
