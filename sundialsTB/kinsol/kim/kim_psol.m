function [ret, flag, new_data] = kim_psol(y, yscale, fy, fscale, v, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%
if isempty(data)
  [ret, flag] = feval(fct,y,yscale,fy,fscale,v);
  new_data = [];
else
  [ret, flag, new_data] = feval(fct,y,yscale,fy,fscale,v,data);
end
