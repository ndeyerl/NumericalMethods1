%with shifts
m=50;
a=[2,-1,zeros(1,m-2)];
A0=toeplitz(a);
A=A0;
n=100;
dVec=[];
sVec=[];
for k=1:n
if k==1
mu=0;
[Q,R]=qr(A-mu*eye(m));
A=R*Q+mu*eye(m);
dVec(:,k)=diag(A);
sVec(:,k)=diag(A,-1);
else
mu=A(m,m);
[Q,R]=qr(A-mu*eye(m));
A=R*Q+mu*eye(m);
dVec(:,k)=diag(A);
sVec(:,k)=diag(A,-1);
end
end
%without shifts
m=50;
a=[2,-1,zeros(1,m-2)];
A0=toeplitz(a);
A=A0;
n=100;
dVec=[];
sVec=[];
for k=1:n
mu=0;
[Q,R]=qr(A-mu*eye(m));
A=R*Q+mu*eye(m);
dVec(:,k)=diag(A);
sVec(:,k)=diag(A,-1);
end
