function [ret, flag, new_data] = cvm_jtv(type, varargin)

%
% Wrapper around the actual user-provided Matlab function
%

switch type
  
 case 1

  % Forward ODE
  
  t    = varargin{1};
  y    = varargin{2};
  fy   = varargin{3};
  v    = varargin{4};
  fct  = varargin{5};
  data = varargin{6};

  if isempty(data)
    [ret, flag] = feval(fct,t,y,fy,v);
    new_data = [];
  else
    [ret, flag, new_data] = feval(fct,t,y,fy,v,data);
  end
  
 case -1
  
  % Backward ODE
  
  t    = varargin{1};
  y    = varargin{2};
  yB   = varargin{3};
  fyB  = varargin{4};
  vB   = varargin{5};
  fct  = varargin{6};
  data = varargin{7};
  
  if isempty(data)
    [ret, flag] = feval(fct,t,y,yB,fyB,vB);
    new_data =[];
  else
    [ret, flag, new_data] = feval(fct,t,y,yB,fyB,vB,data);
  end

  
end