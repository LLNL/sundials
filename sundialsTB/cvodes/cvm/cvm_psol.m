function [z, flag, new_data] = cvm_psol(t, y, fy, r, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [z, flag] = feval(fct,t,y,fy,r);
  new_data = [];
else
  [z, flag, new_data] = feval(fct,t,y,fy,r,data);
end
