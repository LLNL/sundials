function [J, flag, new_data] = idm_bjac(t, yy, yp, rr, cj, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [J, flag] = feval(fct,t,yy,yp,rr,cj);
  new_data = [];
else
  [J, flag, new_data] = feval(fct,t,yy,yp,rr,cj,data);
end

