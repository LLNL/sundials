function [rrS, flag, new_data] = idm_resS(t, yy, yp, rr, Ns, yyS, ypS, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

N = length(yy);
yyS = reshape(yyS, N, Ns);
ypS = reshape(ypS, N, Ns);

if isempty(data)
  [rrS, flag] = feval(fct,t,yy,yp,rr,yyS,ypS);
  new_data = [];
else
  [rrS, flag, new_data] = feval(fct,t,yy,yp,rr,yyS,ypS,data);
end

rrS = reshape(rrS, N*Ns, 1);