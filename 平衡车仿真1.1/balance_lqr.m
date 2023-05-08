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
     0, 0, 0,  0, 100, 0;
     0, 0, 0,  0,   0, 0];

R = [1 0;
     0 1];

K = lqr(A,B,Q,R)


%%
%三维建模real

clc;clear;

syms F ddth ddx th TR TL T0 ddbe dbe be;

g = 9.81;

m = 0.98;
M = 5.55;

r = 0.1;

d = 0.23;%高
l = 0.6;%长

Jy = 1/3*M*d^2; %0.064;
Jz = (l/2)^2*(M + m*2)/2;

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

%%

clc;clear;

syms F ddth ddx th TR TL T0 ddbe dbe be;

g = 9.81;

m = 0.98;
M = 5.55;

r = 0.1;

d = 0.23;%高
l = 0.6;%长

Jy = 1/3*M*d^2; %0.064;
Jz = (l/2)^2*(M + m*2)/2;

I = 1/2*m*r^2;

A1 = -(M*d*r*g)^2/(2*I*Jy+2*I*M*d^2+M*r^2*Jy+2*Jy*m*r^2+2*M*m*d^2*r^2);







