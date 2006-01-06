function [jcur, status, data] = cvkx_pset(t,u,fu,jok,gm,data)
%CVKX_PSET - Preconditioner setup function for the CVKX example problem.
%
%   See also: cvkx, CVPrecSetupFn

% Radu Serban <radu@llnl.gov>
% Copyright (c) 2005, The Regents of the University of California.
% $Revision: 1.1 $Date$

persistent Jbd

mx = data.mx;
my = data.my;

if jok

  % Copy Jbd to P
  
  P = Jbd;
  jcur = false;
  
else

  % Generate Jbd from scratch and copy to P

  xmin = data.xmin;
  dx = data.dx;

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
    q4 = exp(-a4/s);
  else
    q4 = 0.0;
  end

  u = reshape(u,2,mx*my);
  
  for jy = 1:my
    ydn = ymin + (jy - 1.5)*dy;
    yup = ydn + dy;
    cydn = vdco * exp(0.2*ydn);
    cyup = vdco * exp(0.2*yup);
    diag = -(cydn + cyup + 2.0*hdco);
    i = (jy-1)*mx;
    for jx = 1:mx
      ii = i + jx;
      c1 = u(1,ii);
      c2 = u(2,ii);
      Jbd(1,1,ii) = (-q1*c3 - q2*c2) + diag;
      Jbd(1,2,ii) = -q2*c1 + q4;
      Jbd(2,1,ii) = q1*c3 - q2*c2;
      Jbd(2,2,ii) = (-q2*c1 - q4) + diag;
    end
  end
  
  P = Jbd;
  jcur = true;

end

% Scale by -gamma and add identity
P = - gm*P;
for i = 1:mx*my
  P(:,:,i) = eye(2) + P(:,:,i);
end

status = 0;

data.P = P;
