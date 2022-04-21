function [new_mondata] = idm_monitorB(call, idxB, t, y, yQ, fct, mondata)

%
% Wrapper around the actual user-provided Matlab function
%

N = length(y);

idxB = idxB+1;
if isempty(mondata)
  feval(fct, call, idxB, t, y, yQ);
  new_mondata = [];
else
  new_mondata = feval(fct, call, idxB, t, y, yQ, mondata);
end
