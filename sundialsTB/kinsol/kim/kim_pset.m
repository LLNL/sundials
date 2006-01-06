function [status, new_data] = kim_pset(y, yscale, fy, fscale, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  status = feval(fct,y,yscale,fy,fscale);
  new_data =[];
else
  [status, new_data] = feval(fct,y,yscale,fy,fscale,data);
end
