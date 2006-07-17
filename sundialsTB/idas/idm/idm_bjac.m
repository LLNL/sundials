function [ret, flag, new_data] = idm_bjac(type, varargin)

%
% Wrapper around the actual user-provided Matlab function
%

switch type
  
 case 1

  % Forward ODE
  
  t    = varargin{1};
  yy   = varargin{2};
  yp   = varargin{3};
  rr   = varargin{4};
  cj   = varargin{5};
  fct  = varargin{6};
  data = varargin{7};

  if isempty(data)
    [ret, flag] = feval(fct,t,yy,yp,rr,cj);
    new_data = [];
  else
    [ret, flag, new_data] = feval(fct,t,yy,yp,rr,cj,data);
  end
 
 case -1
  
  % Backward ODE
  
  t    = varargin{1};
  yy   = varargin{2};
  yp   = varargin{3};
  yyB  = varargin{4};
  ypB  = varargin{5};
  rrB  = varargin{6};
  cjB  = varargin{7};
  fct  = varargin{8};
  data = varargin{9};
  
  if isempty(data)
    [ret, flag] = feval(fct,t,yy,yp,yyB,ypB,rrB,cjB);
    new_data = [];
  else
    [ret, flag, new_data] = feval(fct,t,yy,yp,yyB,ypB,rrB,cjB,data);
  end
  
end