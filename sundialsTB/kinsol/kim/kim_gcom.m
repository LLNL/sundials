function [new_data] = kim_gcom(y, f, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  feval(fct,y);
  new_data = [];
else
  [new_data] = feval(fct,y,data);
end

