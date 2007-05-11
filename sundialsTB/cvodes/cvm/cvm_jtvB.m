function [JvB, flag, new_data] = cvm_jtvB(t, y, yB, fyB, vB, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [JvB, flag] = feval(fct,t,y,yB,fyB,vB);
  new_data =[];
else
  [JvB, flag, new_data] = feval(fct,t,y,yB,fyB,vB,data);
end
