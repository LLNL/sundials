function [flag, new_data] = cvm_gcomB(t, y, yB, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  flag = feval(fct,t,y,yB);
  new_data = [];
else
  [flag, new_data] = feval(fct,t,y,yB,data);
end
