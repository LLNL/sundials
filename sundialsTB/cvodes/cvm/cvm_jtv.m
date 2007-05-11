function [Jv, flag, new_data] = cvm_jtv(t, y, fy, v, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [Jv, flag] = feval(fct,t,y,fy,v);
  new_data = [];
else
  [Jv, flag, new_data] = feval(fct,t,y,fy,v,data);
end
