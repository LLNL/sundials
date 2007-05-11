function [flB, flag, new_data] = cvm_glocB(t, y, yB, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [flB, flag] = feval(fct,t,y,yB);
  new_data = [];
else
  [flB, flag, new_data] = feval(fct,t,y,yB,data);
end
