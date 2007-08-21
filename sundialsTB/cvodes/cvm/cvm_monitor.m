function new_mondata = cvm_monitor(call, t, y, yQ, Ns, yS, fct, mondata)

%
% Wrapper around the actual user-provided Matlab function
%

N = length(y);

yS = reshape(yS,N,Ns);

if isempty(mondata)
  feval(fct, call, t, y, yQ, yS);
  new_mondata = [];
else
  new_mondata = feval(fct, call, t, y, yQ, yS, mondata);
end

