function [resB, flag, new_data] = idm_glocB(t, yy, yp, yyB, ypB, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [resB, flag] = feval(fct,t,yy,yp,yyB,ypB);
  new_data = [];
else
  [resB, flag, new_data] = feval(fct,t,yy,yp,yyB,ypB,data);
end

