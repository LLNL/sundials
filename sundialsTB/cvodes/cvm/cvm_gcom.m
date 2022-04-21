function [flag, new_data] = cvm_gcom(t, y, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  flag = feval(fct,t,y);
  new_data = [];
else
  [flag, new_data] = feval(fct,t,y,data);
end
