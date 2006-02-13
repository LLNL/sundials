function [ret, flag, new_data] = cvm_djac(type, varargin)

%
% Wrapper around the actual user-provided Matlab function
%

switch type
  
 case 1

  % Forward ODE

  t    = varargin{1};
  y    = varargin{2};
  fy   = varargin{3};
  fct  = varargin{4};
  data = varargin{5};

  if isempty(data)
    [ret, flag] = feval(fct,t,y,fy);
    new_data = [];
  else
    [ret, flag, new_data] = feval(fct,t,y,fy,data);
  end

 case -1
  
  % Backward ODE
  
  t    = varargin{1};
  y    = varargin{2};
  yB   = varargin{3};
  fyB  = varargin{4};
  fct  = varargin{5};
  data = varargin{6};
  
  if isempty(data)
    [ret, flag] = feval(fct,t,y,yB,fyB);
    new_data = [];
  else
    [ret, flag, new_data] = feval(fct,t,y,yB,fyB,data);
  end
  
end