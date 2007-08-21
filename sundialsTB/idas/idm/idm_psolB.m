function [zB, flag, new_data] = idm_psolB(t, yy, yp, yyB, ypB, rrB, rB, cjB, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [zB, flag] = feval(fct,t,yy,yp,yyB,ypB,rrB,rB,cjB);
  new_data = [];
else
  [zB, flag, new_data] = feval(fct,t,yy,yp,yyB,ypB,rrB,rB,cjB,data);
end

