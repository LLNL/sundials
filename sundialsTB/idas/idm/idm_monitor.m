function [new_mondata] = idm_monitor(call, t, yy, yp, yQ, Ns, yyS, ypS, fct, mondata)

%
% Wrapper around the actual user-provided Matlab function
%

N = length(yy);
yyS = reshape(yyS, N, Ns);
ypS = reshape(ypS, N, Ns);

if isempty(mondata)
  feval(fct, call, t, yy, yp, yQ, yyS, ypS);
  new_mondata = [];
else
  new_mondata = feval(fct, call, t, yy, yp, yQ, yyS, ypS, mondata);
end

