function [JvB, flag, new_data] = idm_jtvB(t, yy, yp, yyB, ypB, rrB, vB, cjB, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [JvB, flag] = feval(fct,t,yy,yp,yyB,ypB,rrB,vB,cjB);
  new_data = [];
else
  [JvB, flag, new_data] = feval(fct,t,yy,yp,yyB,ypB,rrB,vB,cjB,data);
end
