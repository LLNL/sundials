function [ret, flag, new_data] = idm_psol(type, varargin)

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
  r    = varargin{5};
  cj   = varargin{6};
  fct  = varargin{7};
  data = varargin{8};

  if isempty(data)
    [ret, flag] = feval(fct,t,yy,yp,rr,r,cj);
    new_data = [];
  else
    [ret, flag, new_data] = feval(fct,t,yy,yp,rr,r,cj,data);
  end
  
 case -1
  
  % Backward ODE
  
  t    = varargin{1};
  yy   = varargin{2};
  yp   = varargin{3};
  yyB  = varargin{4};
  ypB  = varargin{5};
  rrB  = varargin{6};
  rB   = varargin{7};
  cjB  = varargin{8};
  fct  = varargin{9};
  data = varargin{10};
  
  if isempty(data)
    [ret, flag] = feval(fct,t,yy,yp,yyB,ypB,rrB,rB,cjB);
    new_data = [];
  else
    [ret, flag, new_data] = feval(fct,t,yy,yp,yyB,ypB,rrB,rB,cjB,data);
  end
  
end