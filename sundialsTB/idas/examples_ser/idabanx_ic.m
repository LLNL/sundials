function [t,yy,yp,id,cnstr] = idabanx_ic(data)

m = data.m;
N = data.N;
dx = data.dx;

id = ones(N,1);
cnstr = ones(N,1);
yy = zeros(N,1);
yp = zeros(N,1);

t = 0.0;

% Initialize yy on all grid points. */ 
for j=0:m-1
  yfact = dx * j;
  offset = m*j;
  for i=0:m-1
    xfact = dx * i;
    loc = offset + i + 1;
    yy(loc) = 16.0 * xfact * (1.0 - xfact) * yfact * (1.0 - yfact);
  end
end
  
% The residual gives the negative of ODE RHS values at 
% interior points.
yp = zeros(N,1);
[yp,flag,new_data] = idabanx_f(t,yy,yp,data);
yp = -yp;

% Finally, set values of yy, yp, and id at boundary points.
for j=0:m-1
  offset = m*j;
  for i=0:m-1
    loc = offset + i + 1;
    if (j == 0 || j == m-1 || i == 0 || i == m-1 )
        yy(loc) = 0.1;
        yp(loc) = 0.0;
        id(loc) = 0.0;
    end
  end
end
