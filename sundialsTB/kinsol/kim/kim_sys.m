function [ret, flag, new_data] = kim_sys(y, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [ret, flag] = feval(fct,y);
  new_data = [];
else
  [ret, flag, new_data] = feval(fct,y,data);
end

