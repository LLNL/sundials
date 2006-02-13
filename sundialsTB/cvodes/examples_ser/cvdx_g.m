function [g, flag, new_data] = cvdx_g(t,y,data)
%CVDX_G - Root-finding function for the CVDX example problem.
%
%   See also: cvdx, CVRootFn 

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/01/06 18:59:49 $


g(1) = y(1) - 0.0001;
g(2) = y(3) - 0.01;

flag = 0;
new_data = [];
