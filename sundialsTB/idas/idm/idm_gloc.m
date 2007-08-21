function [res, flag, new_data] = idm_gloc(t, yy, yp, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [res, flag] = feval(fct,t,yy,yp);
  new_data = [];
else
  [res, flag, new_data] = feval(fct,t,yy,yp,data);
end

