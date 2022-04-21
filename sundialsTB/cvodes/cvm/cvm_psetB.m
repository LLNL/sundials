function [jcurB, flag, new_data] = cvm_psetB(t, y, yB, fyB, jokB, gmB, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [jcurB, flag] = feval(fct,t,y,yB,fyB,jokB,gmB);
  new_data = [];
else
  [jcurB, flag, new_data] = feval(fct,t,y,yB,fyB,jokB,gmB,data);
end
