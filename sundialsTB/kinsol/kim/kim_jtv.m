function [ret, new_y, flag, new_data] = kim_jtv(y, v, new_y, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [ret, new_y, flag] = feval(fct,y,v,new_y);
  new_data = [];
else
  [ret, new_y, flag, new_data] = feval(fct,y,v,new_y,data);
end
