function [ret, new_data] = kim_sys(y, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  ret = feval(fct,y);
  new_data = [];
else
  [ret, new_data] = feval(fct,y,data);
end
