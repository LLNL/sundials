function [ud, flag, new_data] = cvkx_f(t, u, data)
%CVKX_F - RHS function for the CVKX and CVKXB example problems.
%
%   See also: cvkx, cvkxb, CVRhsFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.2 $Date: 2006/01/06 18:59:49 $


mx = data.mx;
xmin = data.xmin;
dx = data.dx;

my = data.my;
ymin = data.ymin;
dy = data.dy;

om = data.om;
q1 = data.q1;
q2 = data.q2;
c3 = data.c3;
a3 = data.a3;
a4 = data.a4;
hdco = data.hdco;
haco = data.haco;
vdco = data.vdco;

s = sin(om*t);
if s > 0.0
  q3 = exp(-a3/s);
  q4 = exp(-a4/s);
else
  q3 = 0.0;
  q4 = 0.0;
end

u = reshape(u, 2,mx*my);

for jy = 1:my
  ydn = ymin + (jy - 1.5)*dy;
  yup = ydn + dy;
  cydn = vdco * exp(0.2*ydn);
  cyup = vdco * exp(0.2*yup);
  i = (jy-1)*mx;
  idn = -mx;
  if jy == 1
    idn = mx;
  end
  iup = mx;
  if jy == my
    iup = -mx;
  end
  for jx = 1:mx
    ii = i + jx;
    c1 = u(1,ii);
    c2 = u(2,ii);
    % kinetic rate terms
    qq1 = q1 * c1 * c3;
    qq2 = q2 * c1 * c2;
    qq3 = q3 * c3;
    qq4 = q4 *c2;
    rkin1 = -qq1 - qq2 + 2.0*qq3 + qq4;
    rkin2 = qq1 - qq2 - qq4;
    % vertical diffusion
    c1dn = u(1,ii+idn);
    c2dn = u(2,ii+idn);
    c1up = u(1,ii+iup);
    c2up = u(2,ii+iup);
    vertd1 = cyup*(c1up - c1) - cydn*(c1 - c1dn);
    vertd2 = cyup*(c2up - c2) - cydn*(c2 - c2dn);
    % horizontal diffusion and advection
    ileft = -1;
    if jx == 1
      ileft = 1;
    end
    iright = 1;
    if jx == mx
      iright = -1;
    end
    c1lt = u(1,ii+ileft);
    c2lt = u(2,ii+ileft);
    c1rt = u(1,ii+iright);
    c2rt = u(2,ii+iright);
    hord1 = hdco * (c1rt-2.0*c1+c1lt);
    hord2 = hdco * (c2rt-2.0*c2+c2lt);
    horad1 = haco * (c1rt-c1lt);
    horad2 = haco * (c2rt-c2lt);
    % load into ud
    ud(1,ii) = vertd1 + hord1 + horad1 + rkin1; 
    ud(2,ii) = vertd2 + hord2 + horad2 + rkin2;
  end
  
end

ud = reshape(ud,2*mx*my,1);

flag = 0;
new_data = [];

new_data = data;