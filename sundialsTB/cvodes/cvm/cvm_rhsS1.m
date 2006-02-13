function [ySd, flag, new_data] = cvm_rhsS1(is, t, y, yd, yS, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [ySd, flag] = feval(fct,is,t,y,yd,yS);
  new_data = [];
else
  [ySd, flag, new_data] = feval(fct,is,t,y,yd,yS,data);
end
