function [ret, new_y, status, new_data] = kim_jtv(y, v, new_y, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [ret, new_y, status] = feval(fct,y,v,new_y);
  new_data = [];
else
  [ret, new_y, status, new_data] = feval(fct,y,v,new_y,data);
end
