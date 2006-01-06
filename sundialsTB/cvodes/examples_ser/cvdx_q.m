function [qd, new_data] = cvdx_q(t, y, data)
%CVDX_Q - quadrature function for the CVADX example problem.
%
%   See also: cvadx, CVQuadRhsFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$


qd = y(3);

new_data = [];
