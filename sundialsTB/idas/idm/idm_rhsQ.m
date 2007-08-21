function [qd, flag, new_data] = idm_rhsQ(t, yy, yp, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [qd, flag] = feval(fct,t,yy,yp);
  new_data =[];
else
  [qd, flag, new_data] = feval(fct,t,yy,yp,data);
end
