function [flag, new_data] = idm_pset(t, yy, yp, rr, cj, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  flag = feval(fct,t,yy,yp,rr,cj);
  new_data = [];
else
  [flag, new_data] = feval(fct,t,yy,yp,rr,cj,data);
end
