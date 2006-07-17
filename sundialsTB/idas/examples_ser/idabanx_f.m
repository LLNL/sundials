function [rr,flag,new_data] = idabanx_f(t,yy,yp,data)

m = data.m;
N = data.N;
dx = data.dx;
c = data.c;

% Initialize resval to uu, to take care of boundary equations.
rr = yy;
  
% Loop over interior points; set rr = yp - (central difference).
for j = 1:m-2
  offset = m*j;
  for i = 1:m-2
    loc = offset + i + 1;
    rr(loc) = yp(loc) - c * ...
	  (yy(loc-1) + yy(loc+1) + yy(loc-m) + yy(loc+m) - 4.0*yy(loc));
  end
end

flag = 0;
new_data = [];