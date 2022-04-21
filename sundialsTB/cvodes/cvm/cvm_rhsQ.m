function [qd, flag, new_data] = cvm_rhsQ(t, y, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [qd, flag] = feval(fct,t,y);
  new_data =[];
else
  [qd, flag, new_data] = feval(fct,t,y,data);
end

