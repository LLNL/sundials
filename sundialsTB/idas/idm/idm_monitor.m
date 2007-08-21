function new_mondata = idm_monitor(call, t, yy, yQ, Ns, yyS, fct, mondata)

%
% Wrapper around the actual user-provided Matlab function
%

N = length(yy);
yyS = reshape(yyS, N, Ns);

if isempty(mondata)
  feval(fct, call, t, yy, yQ, yyS);
  new_mondata = [];
else
  new_mondata = feval(fct, call, t, yy, yQ, yyS, mondata);
end

