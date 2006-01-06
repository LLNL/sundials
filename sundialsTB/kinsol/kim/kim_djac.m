function [ret, new_data] = kim_djac(y, fy, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  ret = feval(fct,y,fy);
  new_data = [];
else
  [ret, new_data] = feval(fct,y,fy,data);
end
