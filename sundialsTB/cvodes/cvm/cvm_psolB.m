function [zB, flag, new_data] = cvm_psolB(t, y, yB, fyB, rB, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [zB, flag] = feval(fct,t,y,yB,fyB,rB);
  new_data = [];
else
  [zB, flag, new_data] = feval(fct,t,y,yB,fyB,rB,data);
end
