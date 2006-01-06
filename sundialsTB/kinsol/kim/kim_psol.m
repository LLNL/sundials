function [ret, status, new_data] = kim_psol(y, yscale, fy, fscale, v, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%
if isempty(data)
  [ret, status] = feval(fct,y,yscale,fy,fscale,v);
  new_data = [];
else
  [ret, status, new_data] = feval(fct,y,yscale,fy,fscale,v,data);
end
