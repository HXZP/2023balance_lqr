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
d = 0.2;
l = 0.4;
Jy = 1/3*M*d^2; %0.064;
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
     0, 0, 0,  0,   1, 0;
     0, 0, 0,  0,   0, 1000];

R = [1 0;
     0 1];

K = lqr(A,B,Q,R)



%%
%三维建模real

clc;clear;

syms F ddth ddx th TR TL T0 ddbe dbe be T;
% syms g M d l Jy Jz Jx m r I;

g = 9.7833;

M = 6.754+0.88;
d = 0.05244;
l = 0.48;

Jy = 1206539.81/1000/10000;
Jz = 0.35;
Jx = 0.46;

m = 1.423-0.44;
r = 0.1;
I = 30022.65/1000/10000;

f1 = (M*d^2+Jy)*ddth == M*g*d*th-M*d*ddx-F*r;

f2 = F*r^2-M*r^2*d*ddth == (2*I+2*m*r^2+M*r^2)*ddx;

f3 = ddbe == T0/(r*(2*Jz/l+l*(m*r^2+I)/r^2));

[s_ddth,s_ddx,s_ddbe] = solve(f1,f2,f3,ddth,ddx,ddbe)

a = diff(s_ddth,th);%2
b = diff(s_ddth,F);
c = diff(s_ddx,th);%1
e = diff(s_ddx,F);%
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

Q = [500000, 0, 0, 0, 0, 0;
     0, 900000, 0, 0, 0, 0;
     0, 0, 5000000, 0, 0, 0;
     0, 0, 0,  100000, 0, 0;
     0, 0, 0,  0,   500000, 0;
     0, 0, 0,  0,   0, 50000];

R = [500 0;
     0 500];

K = lqr(A,B,Q,R)

Acl = A-B*K

% aaa = -(M*r*r*d*M*g*d)/(M*d*r*r*M*d-(2*I+M*r*r+2*m*r*r)*(M*d*d+Jy))%自己之前算的
% bbb = M*d*r*M*d*g/(r*(M+2*m+2*I/(r*r))*(Jy+M*d*d)-M*r*d*M*d)%其他学校
% ccc = -(M^2*d^2*g*r^2)/(2*I*Jy + 2*I*M*d^2 + Jy*M*r^2 + 2*Jy*m*r^2 + 2*M*d^2*m*r^2)%代码
% ddd = M*M*r*r*d*d*g/((2*I+2*m*r*r+M*r*r)*(M*d*d+Jy)-M*M*r*r*d*d)
%%
%三维建模real

clc;clear;

syms F ddth ddx th TR TL T0 ddbe dbe be T;

% syms g M d l Jy Jz Jx m r I;
% syms a b c e f h k


% g = 9.7833;
% 
% M = 6.754+0.88;
% d = 0.05244;
% l = 0.48;
% 
% Jy = 1206539.81/1000/10000;
% Jz = 0.35;
% Jx = 0.46;
% 
% m = 1.423-0.44;
% r = 0.1;
% I = 30022.65/1000/10000;


g = 9.81;

m = 0.98;
M = 5.55;

r = 0.1;

d = 0.23;%高
l = 0.6;%长

Jy = 1/3*M*d^2; %0.064;
Jz = (l/2)^2*(M + m*2)/2;

I = 1/2*m*r^2;

f1 = (TR+TL)*r-M*r^2*d*ddth == (2*I+2*m*r^2+M*r^2)*ddx;

f2 = (M*d^2+Jy)*ddth == M*g*d*th-M*d*ddx-(TR+TL);

f3 = ddbe == (TR-TL)/(r*(2*Jz/l+l*(m*r^2+I)/r^2));

[s_ddx,s_ddth,s_ddbe] = solve(f1,f2,f3,ddx,ddth,ddbe)



A1 = diff(s_ddx,th);
B1 = diff(s_ddx,TR);
A2 = diff(s_ddth,th);
B2 = diff(s_ddth,TR);
B3 = diff(s_ddbe,TR);

A1 = double(A1);
B1 = double(B1);
A2 = double(A2);
B2 = double(B2);
B3 = double(B3);


A = [0 1 0 0 0 0;
     0 0 A1 0 0 0;
     0 0 0 1 0 0;
     0 0 A2 0 0 0;
     0 0 0 0 0 1;
     0 0 0 0 0 0];

B = [0   0;
     B1  B1;
     0   0;
     B2  B2;
     0   0;
     B3 -B3];

eig(A)

C0 = [B A*B A^2*B A^3*B];

rank(C0)

Q = [1, 0, 0, 0, 0, 0;
     0, 100, 0, 0, 0, 0;
     0, 0, 10, 0, 0, 0;
     0, 0, 0,  100, 0, 0;
     0, 0, 0,  0,   10, 0;
     0, 0, 0,  0,   0, 10];

R = [100 0;
     0 100];

K = lqr(A,B,Q,R)

Acl = A-B*K;

K11 = K(1,1);
K12 = K(1,2);
K13 = K(1,3);
K14 = K(1,4);
K15 = K(1,5);
K16 = K(1,6);

K21 = K(2,1);
K22 = K(2,2);
K23 = K(2,3);
K24 = K(2,4);
K25 = K(2,5);
K26 = K(2,6);


v = 4;

