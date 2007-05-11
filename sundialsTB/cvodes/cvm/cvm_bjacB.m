function [JB, flag, new_data] = cvm_bjacB(t, y, yB, fyB, fct, data)

if isempty(data)
  [JB, flag] = feval(fct,t,y,yB,fyB);
  new_data = [];
else
  [JB, flag, new_data] = feval(fct,t,y,yB,fyB,data);
end
