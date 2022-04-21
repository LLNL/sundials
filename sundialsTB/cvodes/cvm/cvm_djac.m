function [J, flag, new_data] = cvm_djac(t, y, fy, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [J, flag] = feval(fct,t,y,fy);
  new_data = [];
else
  [J, flag, new_data] = feval(fct,t,y,fy,data);
end

