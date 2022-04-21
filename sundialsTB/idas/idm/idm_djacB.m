function [JB, flag, new_data] = idm_djacB(t, yy, yp, yyB, ypB, rrB, cjB, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [JB, flag] = feval(fct,t,yy,yp,yyB,ypB,rrB,cjB);
  new_data = [];
else
  [JB, flag, new_data] = feval(fct,t,yy,yp,yyB,ypB,rrB,cjB,data);
end
