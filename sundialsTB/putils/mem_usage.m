function mem = mem_usage(where)
%
% where can be either 'self' or 'base'
%

if strcmp(where,'self')
  WS = 'caller';
else
  WS = 'base';
end

s = evalin(WS,'whos');
mem = 0;
for i = 1:length(s)
  mem = mem + s(i).bytes;
end

