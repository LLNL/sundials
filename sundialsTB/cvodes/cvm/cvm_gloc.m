function [fl, flag, new_data] = cvm_gloc(t, y, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [fl, flag] = feval(fct,t,y);
  new_data = [];
else
  [fl, flag, new_data] = feval(fct,t,y,data);
end
