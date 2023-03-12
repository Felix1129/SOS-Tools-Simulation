clc
clear
syms x1 x2 v1 v2 ;
vars=[x1;x2;v1;v2];
x=[x1;x2];
v=[v1;v2];
x_h=x;
T=eye(2);
I=eye(2);
prog=sosprogram(vars);

A1=[-1+x1+x1^2+x1*x2-x2^2, 1; -1, -1];
A2=[-1+x1+x1^2+x1*x2-x2^2, 1; 0.2172, -1];
B1=[x1;0];
B2=[x1;0];

[prog,X]=sospolymatrixvar(prog,monomials(x,[0:0]),[2,2],'symmetric');
X
[prog,M1]=sospolymatrixvar(prog,monomials(x,[0:1]),[1,2]);
[prog,M2]=sospolymatrixvar(prog,monomials(x,[0:1]),[1,2]);
M1
[prog,eps1]=sossosvar(prog,monomials(x,[0:0]));
[prog,eps211]=sossosvar(prog,monomials(x,[0:0]));
[prog,eps212]=sossosvar(prog,monomials(x,[0:0]));
[prog,eps222]=sossosvar(prog,monomials(x,[0:0]));


SOS1=v.'*(X-eps1.*I)*v;
prog=sosineq(prog,SOS1,'Mineq');

SOS2=-v.'*(T*A1*X-T*B1*M1+X*A1.'*T.'-M1.'*B1.'*T.'+...
            T*A1*X-T*B1*M1+X*A1.'*T.'-M1.'*B1.'*T.'-...
            diff(X,x2).*(A1(2,:)*x_h)-...
            diff(X,x2).*(A1(2,:)*x_h)+...
            eps211.*I)*v;
prog=sosineq(prog,SOS2,'Mineq');

SOS3=-v.'*(T*A1*X-T*B1*M2+X*A1.'*T.'-M2.'*B1.'*T.'+...
            T*A2*X-T*B2*M1+X*A2.'*T.'-M1.'*B2.'*T.'-...
            diff(X,x2).*(A1(2,:)*x_h)-...
            diff(X,x2).*(A1(2,:)*x_h)+...
            eps212.*I)*v;
prog=sosineq(prog,SOS3,'Mineq');

SOS4=-v.'*(T*A2*X-T*B2*M2+X*A2.'*T.'-M2.'*B2.'*T.'+...
            T*A2*X-T*B2*M2+X*A2.'*T.'-M2.'*B2.'*T.'-...
            diff(X,x2).*(A2(2,:)*x_h)-...
            diff(X,x2).*(A2(2,:)*x_h)+...
            eps212.*I)*v;
prog=sosineq(prog,SOS4,'Mineq');

SOS5=eps1;
prog=sosineq(prog,SOS5,'Mineq');

SOS6=eps211;
prog=sosineq(prog,SOS6,'Mineq');

SOS7=eps212;
prog=sosineq(prog,SOS7,'Mineq');

SOS8=eps222;
prog=sosineq(prog,SOS8,'Mineq');

[prog,info]=sossolve(prog);

X=sosgetsol(prog,X)
M1=sosgetsol(prog,M1)
M2=sosgetsol(prog,M2)
eps1=sosgetsol(prog,eps1)
eps211=sosgetsol(prog,eps211)
eps212=sosgetsol(prog,eps212)
eps222=sosgetsol(prog,eps222)

F1=sosgetsol(prog,M1,15)*inv(sosgetsol(prog,X,15))
F2=sosgetsol(prog,M2,15)*inv(sosgetsol(prog,X,15))
X
P=inv(X)
V=x_h.'*P*x_h
eig_P = eig(P)
u=-(F1*x_h+F2*x_h)
