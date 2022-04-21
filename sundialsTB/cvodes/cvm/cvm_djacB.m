function [JB, flag, new_data] = cvm_djacB(t, y, yB, fyB, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [JB, flag] = feval(fct,t,y,yB,fyB);
  new_data = [];
else
  [JB, flag, new_data] = feval(fct,t,y,yB,fyB,data);
end
