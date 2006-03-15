function [ret, flag, new_data] = kim_bjac(y, fy, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [ret, flag] = feval(fct,y,fy);
  new_data = [];
else
  [ret, flag, new_data] = feval(fct,y,fy,data);
end
