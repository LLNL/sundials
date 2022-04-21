function [flag, new_data] = idm_psetB(t, yy, yp, yyB, ypB, rrB, cjB, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  flag = feval(fct,t,yy,yp,yyB,ypB,rrB,cjB);
  new_data = [];
else
  [flag, new_data] = feval(fct,t,yy,yp,yyB,ypB,rrB,cjBy,data);
end

