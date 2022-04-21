function [Jv, flag, new_data] = idm_jtv(t, yy, yp, rr, v, cj, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [Jv, flag] = feval(fct,t,yy,yp,rr,v,cj);
  new_data = [];
else
  [Jv, flag, new_data] = feval(fct,t,yy,yp,rr,v,cj,data);
end
