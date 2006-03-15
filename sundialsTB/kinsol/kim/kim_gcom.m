function [flag, new_data] = kim_gcom(y, f, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  flag = feval(fct,y);
  new_data = [];
else
  [flag, new_data] = feval(fct,y,data);
end

