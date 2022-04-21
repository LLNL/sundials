function options = idm_options(KeyNames, varargin)

m = length(KeyNames);

% Initialize the output options structure

options = [];
for i = 1:m
  options.(KeyNames{i}) = [];
end

% If the first argument is an options structure, read its non-empty fields
% and update options. Store in j the start of key-value pairs.

arg = varargin{1};

if isa(arg,'struct')
  for i = 1:m
    if isfield(arg,KeyNames{i})
      options.(KeyNames{i}) = arg.(KeyNames{i});
    end
  end
  j = 2;
else
  j = 1;  
end

% The remaining input arguments must be key-value pairs

if rem(nargin-j,2) ~= 0
  error('Arguments must be key-value pairs.');
end

% Process each key-value pair

np = (nargin-j)/2;

keynames = lower(KeyNames);

for i = 1:np
  
  % Get the key
  key = varargin{j};
  
  % key must be a string 
  if ~isstr(key)
    error(sprintf('Argument %d is not a string property name.', j));
  end
  
  % Get the index in keynames that exactly matches the current key
  % (modulo the case)
  ik = strmatch(lower(key), keynames, 'exact');
  if isempty(ik)
    error(sprintf('Unrecognized property "%s"', key));
  end

  % Get the value
  val = varargin{j+1};

  % Set the proper field in options
  options.(KeyNames{ik}) = val;
  
  % move to next pair  
  j = j+2;
  
end
