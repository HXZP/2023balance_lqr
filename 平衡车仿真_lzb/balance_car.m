
clc;clear;

%normal
syms m M r d l g I Jx Jy Jz T;
% m = 1;
% M = 1;
% r = 1;
% d = 1;
% l = 1;
% g = 1;
% I = 1;
% Jx = 1;
% Jy = 1;
% Jz = 1;
% T = 1;

%middle
syms Nr Nl Pr Pl;

%wheel
syms ddxr dxr xr ddar dar ar;
syms ddxl dxl xl ddal dal al;
syms Fr Tr;
syms Fl Tl;
syms ddx;

%pendulum
syms ddth dth th;
syms dvx vx dvy vy;
syms Tp Tn Tg;

% yaw
syms ddbe dbe be;


%1
f1 = Nr == (Tr - I*ddar)/r - ddar*r*m;
f2 = Nl == (Tl - I*ddal)/r - ddal*r*m;

f3 = ddx == (ddar+ddal)*r/2;


%2
f4 = Tp == (Pr + Pl)*d*th;
f5 = Tn == (Nr + Nl)*d;
% f6 = Tg == M*g*d*th;

f7 = ddth == (Tr + Tl - Tp + Tn)/Jz;

f8 = dvx == (Nr+Nl)/M;
f9 = dvy == (Pr+Pl - M*g)/M;

% f10 = ddxr == ddar*r;
% f11 = ddxl == ddal*r;

f12 = dvx == ddx - ddth*d + dth*dth*d*th;
f13 = dvy == ddth*d*th + dth*dth*d;

f14 = Nr*0.5*l - Nl*0.5*l == Jy*ddbe;
f15 = ddbe == (ddar-ddal)*r/l;

% [~,~, a_ddx, ~,~, a_ddth, ~,~, ~,~, ~,a_ddbe] = solve(f1,f2,f3,f4,f5,f7,f8,f9,f12,f13,f14,f15,ddar,ddal,ddx,Tp,Tn,ddth,dvx,dvy,th,dth,ddbe,ddbe)

% [a_ddar,a_ddal,a_ddx] = solve(f1,f2,f3,ddar,ddal,ddx)

% [~,~,~,b_ddth] = solve(f4,f5,f6,f7,Tp,Tn,Tg,ddth)
% 
% [~,~,~,~,a_ddx,a_ddth] = solve(f1,f2,f8,f9,f12,f3,Nr,Nl,dvx,dvy,ddth,ddx)



% fa = l*ddbe/r + ddbe == (Tr-Tl - (Nr-Nl)*r)/(m*r*r+I) + (l/(2*Jy))*(Nr - Nl);
% fb = ddbe == (l/(2*Jy))*(Nr - Nl);
% s_ddbe = solve(fa,ddbe)
% [s_ddbe,~] = solve(fa,fb,ddbe,Jy)


% fc = (Tr - Tl)/r - (m*r*r + I)*l*ddbe/(r*r) == 2*Jy*ddbe/l;
% 
% s_ddbe = solve(fc,ddbe)

%%
clc;clear;

syms Nrl Prl Nr_l;

%middle
syms Nr Nl Pr Pl;


%normal
syms m M r d l g I Jx Jy Jz T;

%wheel
syms ddxr dxr xr ddar dar ar;
syms ddxl dxl xl ddal dal al;
syms Fr Tr;
syms Fl Tl;
syms ddx;

%pendulum
syms ddth dth th;
syms dvx vx dvy vy;
syms Tp Tn Tg;

% yaw
syms ddbe dbe be;



%z
f1 = Nrl == M*(ddx - ddth*d + dth*dth*d*th);
f2 = Prl == M*(ddth*d*th + dth*dth*d) + M*g;

f3 = Prl == Tp/(d*th);
f4 = Nrl == Tn/d;
% f5 = M*g*d*th == Tg;
f6 = ddth == (Tr + Tl - Tp + Tn)/Jz;

[~,~,~,~,a_ddth] = solve(f1,f2,f3,f4,f6,Nrl,Prl,Tp,Tn,ddth)

%x
f7 = Nrl == M*(ddx-ddth*d+dth*dth*th);
f8 = ddx == (r/(2*(m*r*r+I)))*(Tr+Tl) - r*r/(2*(m*r*r+I))*(Nrl);

[a_ddx,~] = solve(f7,f8,ddx,Nrl)


%yaw
f9 = Nr == (Tr - I*ddar)/r - ddar*r*m;
f10 = Nl == (Tl - I*ddal)/r - ddal*r*m;

f11 = (Nr_l)*0.5*l == Jy*ddbe;
f12 = ddbe == (ddar-ddal)*r/l; 
f13 = Nr_l == (Tr - I*ddar)/r - (Tl - I*ddal)/r + ddal*r*m - ddar*r*m;

[a_ddbe,~,~] = solve(f11,f12,f13,ddbe,Nr_l,ddal)

% [a_ddbe,~] = solve(f11,f12,f13,ddal,ddar)


%%
%二维建模

clc;clear;

syms F ddth ddx th;

m = 4;
M = 1.6;
r = 0.1;
d = 0.2;
J = 1/3*M*d^2; %0.064;
g = 9.81;
I = 1/2*m*r^2;

f1 = (M*d^2+J)*ddth == M*g*d*th-M*d*ddx-F*r;
f2 = F*r^2-M*r^2*d*ddth == (I+m*r^2+M*r^2)*ddx;


[s_ddth,s_ddx] = solve(f1,f2,ddth,ddx)

a = diff(s_ddth,th);
b = diff(s_ddth,F);
c = diff(s_ddx,th);
e = diff(s_ddx,F);

a = double(a)
b = double(b)
c = double(c)
e = double(e)



A = [0 1 0 0;
     0 0 c 0;
     0 0 0 1;
     0 0 a 0];
B = [0;e;0;b];

eig(A)

C0 = [B A*B A^2*B A^3*B];

rank(C0)

Q = [0.01, 0, 0, 0;
     0, 1000, 0, 0;
     0, 0, 1000, 0;
     0, 0, 0, 100];

R = 1;

K = lqr(A,B,Q,R)

%%
%三维建模

clc;clear;

syms F ddth ddx th TR TL T0 ddbe dbe be;

g = 9.81;

m = 2;
M = 1.6;
r = 0.1;
d = 0.4;
l = 0.4;
Jy = 1/3*M*d^2;
Jz = (l/2)^2*(M+m*2)/2;

I = 1/2*m*r^2;


f1 = (M*d^2+Jy)*ddth == M*g*d*th-M*d*ddx-F*r;

f2 = F*r^2-M*r^2*d*ddth == (2*I+2*m*r^2+M*r^2)*ddx;

f3 = ddbe == T0/(r*(2*Jz/l+l*(m*r^2+I)/r^2));

[s_ddth,s_ddx,s_ddbe] = solve(f1,f2,f3,ddth,ddx,ddbe);



a = diff(s_ddth,th);
b = diff(s_ddth,F);
c = diff(s_ddx,th);
e = diff(s_ddx,F);
i = diff(s_ddbe,T0);

a = double(a);
b = double(b);
c = double(c);
e = double(e);
i = double(i);


A = [0 1 0 0 0 0;
     0 0 c 0 0 0;
     0 0 0 1 0 0;
     0 0 a 0 0 0;
     0 0 0 0 0 1;
     0 0 0 0 0 0];
B = [0 0;
     e 0;
     0 0;
     b 0;
     0 0;
     0 i];

eig(A)

C0 = [B A*B A^2*B A^3*B];

rank(C0)

Q = [0.01, 0, 0, 0, 0, 0;
     0, 1000, 0, 0, 0, 0;
     0, 0, 1000, 0, 0, 0;
     0, 0, 0,  100, 0, 0;
     0, 0, 0,  0, 100, 0;
     0, 0, 0,  0,   0, 0];

R = [1 0;
     0 1];

K = lqr(A,B,Q,R)

v = 1;
yaw = 20/180*pi;