function [flag, new_data] = idm_gcom(t, yy, yp, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  flag = feval(fct,t,yy,yp);
  new_data = [];
else
  [flag, new_data] = feval(fct,t,yy,yp,data);
end
