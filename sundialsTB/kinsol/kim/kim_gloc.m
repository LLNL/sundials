function [gval, new_data] = kim_gloc(y, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  gval = feval(fct,y);
  new_data = [];
else
  [gval, new_data] = feval(fct,y,data);
end
