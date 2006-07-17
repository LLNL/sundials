function [g, flag, new_data] = idm_root(t, yy, yp, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [g, flag] = feval(fct,t,yy,yp);
  new_data = [];
else
  [g, flag, new_data] = feval(fct,t,yy,yp,data);
end

