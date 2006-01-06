function [JB, new_data] = cvdx_JB(t, y, yB, fyB, data)
%CVDX_JB - adjoint Jacobian for the CVADX example problem.
% 
%   See also: cvadx, CVDenseJacFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$


J = cvdx_J(t,y,[],data);
JB = -J';

new_data = [];