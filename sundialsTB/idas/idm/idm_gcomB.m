function [flag, new_data] = idm_gcomB(t, yy, yp, yyB, ypB, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  flag = feval(fct,t,yy,yp,yyB,ypB);
  new_data = [];
else
  [flag, new_data] = feval(fct,t,yy,yp,yyB,ypB,data);
end
