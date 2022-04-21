function [ySd, flag, new_data] = cvm_rhsS(t, y, yd, Ns, yS, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

N = length(y);
yS = reshape(yS,N,Ns);

if isempty(data)
  [ySd, flag] = feval(fct,t,y,yd,yS);
  new_data = [];
else
  [ySd, flag, new_data] = feval(fct,t,y,yd,yS,data);
end

ySd = reshape(ySd,N*Ns,1);