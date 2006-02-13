function [ret, flag, new_data] = cvm_pset(type, varargin)

%
% Wrapper around the actual user-provided Matlab function
%

switch type
  
 case 1

  % Forward ODE
  
  t    = varargin{1};
  y    = varargin{2};
  fy   = varargin{3};
  jok  = varargin{4};
  gm   = varargin{5};
  fct  = varargin{6};
  data = varargin{7};

  if isempty(data)
    [ret, flag] = feval(fct,t,y,fy,jok,gm);
    new_data = [];
  else
    [ret, flag, new_data] = feval(fct,t,y,fy,jok,gm,data);
  end
  
 case -1
  
  % Backward ODE
  
  t    = varargin{1};
  y    = varargin{2};
  yB   = varargin{3};
  fyB  = varargin{4};
  jokB = varargin{5};
  gmB  = varargin{6};
  fct  = varargin{7};
  data = varargin{8};
  
  if isempty(data)
    [ret, flag] = feval(fct,t,y,yB,fyB,jokB,gmB);
    new_data = [];
  else
    [ret, flag, new_data] = feval(fct,t,y,yB,fyB,jokB,gmB,data);
  end
  
end