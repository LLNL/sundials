function mcvsPollut_FSA_dns()
%mcvsPollut_FSA_dns - Air pollution model
%   J.G Verwer - Gauss-Seidel Iteration for Stiff ODEs from Chemical Kinetics

% Radu Serban <radu@llnl.gov>
% SUNDIALS Copyright Start
% Copyright (c) 2002-2020, Lawrence Livermore National Security
% and Southern Methodist University.
% All rights reserved.
%
% See the top-level LICENSE and NOTICE files for details.
%
% SPDX-License-Identifier: BSD-3-Clause
% SUNDIALS Copyright End
% $Revision$Date: 2007/10/26 16:30:48 $


t0 = 0.0;
tf = 1.0;


Ny = 20;
Np = 25;

plist = [21;22;23;24;25];
Ns = length(plist);

% -------------------
% User data structure
% -------------------

data.p  = [0.35    ; 0.266e2  ; 0.123e5  ; 0.86e-3  ; 0.82e-3  ; ...
           0.15e5  ; 0.13e-3  ; 0.24e5   ; 0.165e5  ; 0.9e4    ; ...
           0.22e-1 ; 0.12e5   ; 0.188e1  ; 0.163e5  ; 0.48e7   ; ...
           0.35e-3 ; 0.175e-1 ; 0.1e9    ; 0.444e12 ; 0.124e4  ; ...
           0.21e1  ; 0.578e1  ; 0.474e-1 ; 0.178e4  ; 0.312e1];

data.plist = plist;

% ---------------------
% CVODES initialization
% ---------------------

options = CVodeSetOptions('UserData',data,...
                          'RelTol',1.e-5,...
                          'AbsTol',1.e-8,...
                          'LinearSolver','Dense');

y0 = [0 ; 0.2   ; 0   ; 0.04 ; 0 ; ...
      0 ; 0.1   ; 0.3 ; 0.01 ; 0 ; ...
      0 ; 0     ; 0   ; 0    ; 0 ; ...
      0 ; 0.007 ; 0   ; 0    ; 0];


CVodeInit(@rhsfn, 'BDF', 'Newton', t0, y0, options);


% ------------------
% FSA initialization
% ------------------

yS0 = zeros(Ny,Ns);

FSAoptions = CVodeSensSetOptions('method','Simultaneous',...
                                 'ErrControl', true,...
                                 'ParamScales', data.p(plist));

CVodeSensInit(Ns, @rhsSfn, yS0, FSAoptions);

% ----------------
% Problem solution
% ----------------

time(1,1) = t0;
sol(1,:) = y0';
sens(1,:,:) = yS0;

t = t0;
it = 1;
while t<tf
  it = it+1;
  [status, t, y, yS] = CVode(tf,'OneStep');
%  [status, t, y] = CVode(tf,'OneStep');
  time(it,1) = t;
  sol(it,:) = y';
  sens(it,:,:) = yS;
end

si = CVodeGetStats

% -------------
% Plot solution
% -------------

figure
hold on
for i = 1:Ny
  plot(time, sol(:,i))
end
set(gca,'XLim',[t0 tf]);
xlabel('time')
grid on
box on

% -----------
% Free memory
% -----------

CVodeFree;

% ===========================================================================

function [yd, flag, new_data] = rhsfn(t, y, data)
% Right-hand side function

p  = data.p;

p1  = p(1);  p2  = p(2);  p3  = p(3);  p4  = p(4);  p5  = p(5);
p6  = p(6);  p7  = p(7);  p8  = p(8);  p9  = p(9);  p10 = p(10);
p11 = p(11); p12 = p(12); p13 = p(13); p14 = p(14); p15 = p(15);
p16 = p(16); p17 = p(17); p18 = p(18); p19 = p(19); p20 = p(20);
p21 = p(21); p22 = p(22); p23 = p(23); p24 = p(24); p25 = p(25);

y1  = y(1);  y2  = y(2);  y3  = y(3);  y4  = y(4);  y5  = y(5);
y6  = y(6);  y7  = y(7);  y8  = y(8);  y9  = y(9);  y10 = y(10);
y11 = y(11); y12 = y(12); y13 = y(13); y14 = y(14); y15 = y(15);
y16 = y(16); y17 = y(17); y18 = y(18); y19 = y(19); y20 = y(20);


r1  = p1*y1;
r2  = p2*y2*y4;
r3  = p3*y2*y5;
r4  = p4*y7;
r5  = p5*y7;

r6  = p6*y6*y7;
r7  = p7*y9;
r8  = p8*y6*y9;
r9  = p9*y2*y11;
r10 = p10*y1*y11;

r11 = p11*y13;
r12 = p12*y2*y10;
r13 = p13*y14;
r14 = p14*y1*y6;
r15 = p15*y3;

r16 = p16*y4;
r17 = p17*y4;
r18 = p18*y16;
r19 = p19*y16;
r20 = p20*y6*y17;

r21 = p21*y19;
r22 = p22*y19;
r23 = p23*y1*y4;
r24 = p24*y1*y19;
r25 = p25*y20;


f1  = -r1-r10-r14-r23-r24+r2+r3+r9+r11+r12+r22+r25;
f2  = -r2-r3-r9-r12+r1+r21;
f3  = -r15+r1+r17+r19+r22;
f4  = -r2-r16-r17-r23+r15;
f5  = -r3+2*r4+r6+r7+r13+r20;

f6  = -r6-r8-r14-r20+r3+2*r18;
f7  = -r4-r5-r6+r13;
f8  = r4+r5+r6+r7;
f9  = -r7-r8;
f10 = -r12+r7+r9;

f11 = -r9-r10+r8+r11;
f12 = r9;
f13 = -r11+r10;
f14 = -r13+r12;
f15 = r14;

f16 = -r18-r19+r16;
f17 = -r20;
f18 = r20;
f19 = -r21-r22-r24+r23+r25;
f20 = -r25+r24;

yd = [f1  ; f2  ; f3  ; f4  ; f5  ; ...
      f6  ; f7  ; f8  ; f9  ; f10 ; ...
      f11 ; f12 ; f13 ; f14 ; f15 ; ...
      f16 ; f17 ; f18 ; f19 ; f20];


flag = 0;
new_data = [];

return

% ===========================================================================

function [ySd, flag, new_data] = rhsSfn(t,y,yd,yS,data)
% Sensitivity right-hand side function

J = Jmat(y,data);
K = Kmat(y,data);

pl = data.plist;

ySd = J*yS + K(:,pl);

flag = 0;
new_data = [];

return

% ===========================================================================

function [J, flag, new_data] = djacfn(t, y, fy, data)

J = Jmat(y, data);

flag = 0;
new_data = [];

return

% ===========================================================================

function J = Jmat(y, data)

p  = data.p;

p1  = p(1);  p2  = p(2);  p3  = p(3);  p4  = p(4);  p5  = p(5);
p6  = p(6);  p7  = p(7);  p8  = p(8);  p9  = p(9);  p10 = p(10);
p11 = p(11); p12 = p(12); p13 = p(13); p14 = p(14); p15 = p(15);
p16 = p(16); p17 = p(17); p18 = p(18); p19 = p(19); p20 = p(20);
p21 = p(21); p22 = p(22); p23 = p(23); p24 = p(24); p25 = p(25);

y1  = y(1);  y2  = y(2);  y3  = y(3);  y4  = y(4);  y5  = y(5);
y6  = y(6);  y7  = y(7);  y8  = y(8);  y9  = y(9);  y10 = y(10);
y11 = y(11); y12 = y(12); y13 = y(13); y14 = y(14); y15 = y(15);
y16 = y(16); y17 = y(17); y18 = y(18); y19 = y(19); y20 = y(20);

% Jacobian obtained symbolically :-)

J = zeros(20, 20);

J( 1,  1) = -p1-p10*y11-p14*y6-p23*y4-p24*y19;
J( 1,  2) = p2*y4+p3*y5+p9*y11+p12*y10;
J( 1,  4) = -p23*y1+p2*y2;
J( 1,  5) = p3*y2;
J( 1,  6) = -p14*y1;
J( 1, 10) = p12*y2;
J( 1, 11) = -p10*y1+p9*y2;
J( 1, 13) = p11;
J( 1, 19) = -p24*y1+p22;
J( 1, 20) = p25;
J( 2,  1) = p1;
J( 2,  2) = -p2*y4-p3*y5-p9*y11-p12*y10;
J( 2,  4) = -p2*y2;
J( 2,  5) = -p3*y2;
J( 2, 10) = -p12*y2;
J( 2, 11) = -p9*y2;
J( 2, 19) = p21;
J( 3,  1) = p1;
J( 3,  3) = -p15;
J( 3,  4) = p17;
J( 3, 16) = p19;
J( 3, 19) = p22;
J( 4,  1) = -p23*y4;
J( 4,  2) = -p2*y4;
J( 4,  3) = p15;
J( 4,  4) = -p2*y2-p16-p17-p23*y1;
J( 5,  2) = -p3*y5;
J( 5,  5) = -p3*y2;
J( 5,  6) = p6*y7+p20*y17;
J( 5,  7) = 2*p4+p6*y6;
J( 5,  9) = p7;
J( 5, 14) = p13;
J( 5, 17) = p20*y6;
J( 6,  1) = -p14*y6;
J( 6,  2) = p3*y5;
J( 6,  5) = p3*y2;
J( 6,  6) = -p6*y7-p8*y9-p14*y1-p20*y17;
J( 6,  7) = -p6*y6;
J( 6,  9) = -p8*y6;
J( 6, 16) = 2*p18;
J( 6, 17) = -p20*y6;
J( 7,  6) = -p6*y7;
J( 7,  7) = -p4-p5-p6*y6;
J( 7, 14) = p13;
J( 8,  6) = p6*y7;
J( 8,  7) = p4+p5+p6*y6;
J( 8,  9) = p7;
J( 9,  6) = -p8*y9;
J( 9,  9) = -p7-p8*y6;
J(10,  2) = -p12*y10+p9*y11;
J(10,  9) = p7;
J(10, 10) = -p12*y2;
J(10, 11) = p9*y2;
J(11,  1) = -p10*y11;
J(11,  2) = -p9*y11;
J(11,  6) = p8*y9;
J(11,  9) = p8*y6;
J(11, 11) = -p9*y2-p10*y1;
J(11, 13) = p11;
J(12,  2) = p9*y11;
J(12, 11) = p9*y2;
J(13,  1) = p10*y11;
J(13, 11) = p10*y1;
J(13, 13) = -p11;
J(14,  2) = p12*y10;
J(14, 10) = p12*y2;
J(14, 14) = -p13;
J(15,  1) = p14*y6;
J(15,  6) = p14*y1;
J(16,  4) = p16;
J(16, 16) = -p18-p19;
J(17,  6) = -p20*y17;
J(17, 17) = -p20*y6;
J(18,  6) = p20*y17;
J(18, 17) = p20*y6;
J(19,  1) = -p24*y19+p23*y4;
J(19,  4) = p23*y1;
J(19, 19) = -p21-p22-p24*y1;
J(19, 20) = p25;
J(20,  1) = p24*y19;
J(20, 19) = p24*y1;
J(20, 20) = -p25;

return

% ===========================================================================

function K = Kmat(y, data)

p  = data.p;

p1  = p(1);  p2  = p(2);  p3  = p(3);  p4  = p(4);  p5  = p(5);
p6  = p(6);  p7  = p(7);  p8  = p(8);  p9  = p(9);  p10 = p(10);
p11 = p(11); p12 = p(12); p13 = p(13); p14 = p(14); p15 = p(15);
p16 = p(16); p17 = p(17); p18 = p(18); p19 = p(19); p20 = p(20);
p21 = p(21); p22 = p(22); p23 = p(23); p24 = p(24); p25 = p(25);

y1  = y(1);  y2  = y(2);  y3  = y(3);  y4  = y(4);  y5  = y(5);
y6  = y(6);  y7  = y(7);  y8  = y(8);  y9  = y(9);  y10 = y(10);
y11 = y(11); y12 = y(12); y13 = y(13); y14 = y(14); y15 = y(15);
y16 = y(16); y17 = y(17); y18 = y(18); y19 = y(19); y20 = y(20);

K = zeros(20, 25);

K( 1,  1) = -y1;
K( 1,  2) = y2*y4;
K( 1,  3) = y2*y5;
K( 1,  9) = y2*y11;
K( 1, 10) = -y1*y11;
K( 1, 11) = y13;
K( 1, 12) = y2*y10;
K( 1, 14) = -y1*y6;
K( 1, 22) = y19;
K( 1, 23) = -y1*y4;
K( 1, 24) = -y1*y19;
K( 1, 25) = y20;
K( 2,  1) = y1;
K( 2,  2) = -y2*y4;
K( 2,  3) = -y2*y5;
K( 2,  9) = -y2*y11;
K( 2, 12) = -y2*y10;
K( 2, 21) = y19;
K( 3,  1) = y1;
K( 3, 15) = -y3;
K( 3, 17) = y4;
K( 3, 19) = y16;
K( 3, 22) = y19;
K( 4,  2) = -y2*y4;
K( 4, 15) = y3;
K( 4, 16) = -y4;
K( 4, 17) = -y4;
K( 4, 23) = -y1*y4;
K( 5,  3) = -y2*y5;
K( 5,  4) = 2*y7;
K( 5,  6) = y6*y7;
K( 5,  7) = y9;
K( 5, 13) = y14;
K( 5, 20) = y6*y17;
K( 6,  3) = y2*y5;
K( 6,  6) = -y6*y7;
K( 6,  8) = -y6*y9;
K( 6, 14) = -y1*y6;
K( 6, 18) = 2*y16;
K( 6, 20) = -y6*y17;
K( 7,  4) = -y7;
K( 7,  5) = -y7;
K( 7,  6) = -y6*y7;
K( 7, 13) = y14;
K( 8,  4) = y7;
K( 8,  5) = y7;
K( 8,  6) = y6*y7;
K( 8,  7) = y9;
K( 9,  7) = -y9;
K( 9,  8) = -y6*y9;
K(10,  7) = y9;
K(10,  9) = y2*y11;
K(10, 12) = -y2*y10;
K(11,  8) = y6*y9;
K(11,  9) = -y2*y11;
K(11, 10) = -y1*y11;
K(11, 11) = y13;
K(12,  9) = y2*y11;
K(13, 10) = y1*y11;
K(13, 11) = -y13;
K(14, 12) = y2*y10;
K(14, 13) = -y14;
K(15, 14) = y1*y6;
K(16, 16) = y4;
K(16, 18) = -y16;
K(16, 19) = -y16;
K(17, 20) = -y6*y17;
K(18, 20) = y6*y17;
K(19, 21) = -y19;
K(19, 22) = -y19;
K(19, 23) = y1*y4;
K(19, 24) = -y1*y19;
K(19, 25) = y20;
K(20, 24) = y1*y19;
K(20, 25) = -y20;

return