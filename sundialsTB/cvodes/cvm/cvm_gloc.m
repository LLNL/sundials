function [ret, flag, new_data] = cvm_gloc(type, varargin)

%
% Wrapper around the actual user-provided Matlab function
%

switch type
  
 case 1

  % Forward ODE
  
  t    = varargin{1};
  y    = varargin{2};
  fct  = varargin{3};
  data = varargin{4};
  
  if isempty(data)
    [ret, flag] = feval(fct,t,y);
    new_data = [];
  else
    [ret, flag, new_data] = feval(fct,t,y,data);
  end

 case -1
  
  % Backward ODE
  
  t    = varargin{1};
  y    = varargin{2};
  yB   = varargin{3};
  fct  = varargin{4};
  data = varargin{5};
  
  if isempty(data)
    [ret, flag] = feval(fct,t,y,yB);
    new_data = [];
  else
    [ret, flag, new_data] = feval(fct,t,y,yB,data);
  end

end