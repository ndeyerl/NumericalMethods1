function [u,x,y,eps]=poisson(n,nIter,omega)
h=1/n;
x=linspace(0,1,n+1);
b=[];
for i=1:n
for j=1:n
b(i,j)=(h.^2)*fcn(x(i),x(j));
end
end
u=zeros(n,n);
eps=[];
for k=1:nIter
uold=u;
for i=2:n-1
for j=2:n-1
u(i,j)=(1-omega)*u(i,j)+omega*((1/4)*(b(i,j)+u(i-1,j)+u(i+1,j)
+u(i,j-1)+u(i,j+1)));
u(i,j)=u(i,j)*inGeom(i,j,n);
end
end
eps(k)=max(max(abs(u-uold)));
end
[x,y]=meshgrid(0:h:1-h,0:h:1-h);
function f=fcn(x,y)
if (x>=.2) && (x<=.4) && (y>=.1) && (y<=.5)
f=1;
elseif (x>=.1) && (x<=.9) && (y>=.8) && (y<=.9)
f=1;
else
f=0;
end
end
function m=inGeom(x,y,n)
x=x/n;
y=y/n;
if (x>=.4) && (x<=.8) && (y>=.3) && (y<=.6)
m=0;
else
m=1;
end
end
end%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
%
Part 1.
[u,x,y,eps]=poisson(20,200,1.1);
subplot(3,1,1); mesh(x,y,u)
hold on
[u,x,y,eps]=poisson(40,400,1.1);
subplot(3,1,2); mesh(x,y,u)
hold on
[u,x,y,eps]=poisson(80,800,1.1);
subplot(3,1,3); mesh(x,y,u)
axis([0 1 0 1 0 .01])
Part 2.
counter=0;
eVec=[];
for omega=1:.05:2
counter=counter+1;
[u,x,y,eps]=poisson(40,400,omega);
eps=eps';
eVec(:,counter)=eps;
end
semilogy(eVec)
