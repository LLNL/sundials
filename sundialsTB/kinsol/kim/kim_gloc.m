function [gval, new_data] = kim_gloc(varargin)

%
% Wrapper around the actual user-provided Matlab function
%

y    = varargin{1};
fct  = varargin{2};
data = varargin{3};

if isempty(data)
  gval = feval(fct,y);
  new_data = [];
else
  [gval, new_data] = feval(fct,y,data);
end
