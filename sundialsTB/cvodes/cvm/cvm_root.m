function [g, flag, new_data] = cvm_root(t, y, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [g, flag] = feval(fct,t,y);
  new_data = [];
else
  [g, flag, new_data] = feval(fct,t,y,data);
end

