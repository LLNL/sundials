function [] = cvm_monitor(call, t, y, yQ, yS, fct, mondata)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(mondata)
  feval(fct, call, t, y, yQ, yS);
else
  feval(fct, call, t, y, yQ, yS, mondata);
end

