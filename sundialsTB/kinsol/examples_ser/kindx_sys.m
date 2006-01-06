function fy = kindx_sys(y)

fy(1) = y(1)^2 + y(2)^2 - 1.0;
fy(2) = y(2) - y(1)^2;

