function [yd, flag, new_data] = cvm_rhs(t, y, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [yd, flag] = feval(fct,t,y);
  new_data = [];
else
  [yd, flag, new_data] = feval(fct,t,y,data);
end

