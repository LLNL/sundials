function [ret, flag, new_data] = cvm_psol(type, varargin)

%
% Wrapper around the actual user-provided Matlab function
%

switch type
  
 case 1

  % Forward ODE
  
  t    = varargin{1};
  y    = varargin{2};
  fy   = varargin{3};
  r    = varargin{4};
  fct  = varargin{5};
  data = varargin{6};

  if isempty(data)
    [ret, flag] = feval(fct,t,y,fy,r);
    new_data = [];
  else
    [ret, flag, new_data] = feval(fct,t,y,fy,r,data);
  end
  
 case -1
  
  % Backward ODE
  
  t    = varargin{1};
  y    = varargin{2};
  yB   = varargin{3};
  fyB  = varargin{4};
  rB   = varargin{5};
  fct  = varargin{6};
  data = varargin{7};
  
  if isempty(data)
    [ret, flag] = feval(fct,t,y,yB,fyB,rB);
    new_data = [];
  else
    [ret, flag, new_data] = feval(fct,t,y,yB,fyB,rB,data);
  end
  
end