function [new_data] = kim_gcom(varargin)

%
% Wrapper around the actual user-provided Matlab function
%

y    = varargin{1};
fct  = varargin{2};
data = varargin{3};
  
if isempty(data)
  feval(fct,y);
  new_data = [];
else
  [new_data] = feval(fct,y,data);
end

