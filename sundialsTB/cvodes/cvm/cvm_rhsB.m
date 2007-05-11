function [yBd, flag, new_data] = cvm_rhsB(type, varargin)
 
%
% Wrapper around the actual user-provided Matlab function
%
 
switch type
   
 case 0
 
  % Not dependent on yS
   
  t    = varargin{1};
  y    = varargin{2};
  yB   = varargin{3};
  fct  = varargin{4};
  data = varargin{5};
   
  if isempty(data)
    [yBd, flag] = feval(fct,t,y,yB);
    new_data = [];
  else
    [yBd, flag, new_data] = feval(fct,t,y,yB,data);
  end
 
 case 1
   
  % Dependent on yS
   
  t    = varargin{1};
  y    = varargin{2};
  yS   = varargin{3};
  yB   = varargin{4};
  fct  = varargin{5};
  data = varargin{6};
   
  if isempty(data)
    [yBd, flag] = feval(fct,t,y,yS,yB);
    new_data = [];
  else
    [yBd, flag, new_data] = feval(fct,t,y,yS,yB,data);
  end
 
end
