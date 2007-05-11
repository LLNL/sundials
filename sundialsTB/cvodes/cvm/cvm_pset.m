function [jcur, flag, new_data] = cvm_pset(t, y, fy, jok, gm, fct, data)

%
% Wrapper around the actual user-provided Matlab function
%

if isempty(data)
  [jcur, flag] = feval(fct,t,y,fy,jok,gm);
  new_data = [];
else
  [jcur, flag, new_data] = feval(fct,t,y,fy,jok,gm,data);
end
