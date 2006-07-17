function [ret, flag, new_data] = idm_gloc(type, varargin)

%
% Wrapper around the actual user-provided Matlab function
%

switch type
  
 case 1

  % Forward ODE
  
  t    = varargin{1};
  yy   = varargin{2};
  yp   = varargin{3};
  fct  = varargin{4};
  data = varargin{5};
  
  if isempty(data)
    [ret, flag] = feval(fct,t,yy,yp);
    new_data = [];
  else
    [ret, flag, new_data] = feval(fct,t,yy,yp,data);
  end

 case -1
  
  % Backward ODE
  
  t    = varargin{1};
  yy   = varargin{2};
  yp   = varargin{3};
  yyB  = varargin{4};
  ypB  = varargin{5};
  fct  = varargin{6};
  data = varargin{7};
  
  if isempty(data)
    [ret, flag] = feval(fct,t,yy,yp,yyB,ypB);
    new_data = [];
  else
    [ret, flag, new_data] = feval(fct,t,yy,yp,yyB,ypB,data);
  end

end