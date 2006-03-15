function [fy, flag] = kindx_sys(y)
%KINDX_SYS - system function for the KINDX example problem.
%
%   See also: kindx, KINsysFn


fy(1) = y(1)^2 + y(2)^2 - 1.0;
fy(2) = y(2) - y(1)^2;

flag = 0; % success

