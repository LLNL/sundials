function [qBd, flag, new_data] = cvm_rhsQB(type, varargin)

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
    [qBd, flag] = feval(fct,t,y,yB);
    new_data = [];
  else
    [qBd, flag, new_data] = feval(fct,t,y,yB,data);
  end

 case 1
 
  % Dependent on yS

  t    = varargin{1};
  y    = varargin{2};
  Ns   = varargin{3};
  yS   = varargin{4};
  yB   = varargin{5};
  fct  = varargin{6};
  data = varargin{7};
  
  N = length(y);
  yS = reshape(yS,N,Ns);
  
  if isempty(data)
    [qBd, flag] = feval(fct,t,y,yS,yB);
    new_data = [];
  else
    [qBd, flag, new_data] = feval(fct,t,y,yS,yB,data);
  end

  
end