%with shifts
m=50;
a=[2,-1,zeros(1,m-2)];
A0=toeplitz(a);
A=A0;
dVec=[];
sVec=[];
tol=1e-15; %epsilon machine
counter=0;
while m>1
if counter==0
mu=0; %first iteration, need to do mu=0
[Q,R]=qr(A-mu*eye(m));
A=R*Q+mu*eye(m);
counter=counter+1; %add one to exit this loop
else
if abs(A(m,m-1))>tol
mu=A(m,m);
[Q,R]=qr(A-mu*eye(m));
A=R*Q+mu*eye(m);
counter=counter+1;
dVec(1:m,counter)=diag(A);
sVec(1:m-1,counter)=diag(A,-1);
else %the A(m,m-1) entry is below tolerance here
counter=counter+1; %add one to make space in dVec/sVec
dVec(1:m,counter)=diag(A);
sVec(1:m-1,counter)=diag(A,-1);
m=m-1;
A=A(1:m,1:m);
end
end
end
%without shifts
m=50;
a=[2,-1,zeros(1,m-2)];
A0=toeplitz(a);
A=A0;
dVec=[];
sVec=[];
tol=1e-15; %epsilon machine
counter=0;
while m>1
if abs(A(m,m-1))>tol
mu=0;
[Q,R]=qr(A-mu*eye(m));
A=R*Q+mu*eye(m);
counter=counter+1;
dVec(1:m,counter)=diag(A);
sVec(1:m-1,counter)=diag(A,-1);
else %the A(m,m-1) entry is below tolerance here
counter=counter+1; %add one to make space in dVec/sVec
dVec(1:m,counter)=diag(A);
sVec(1:m-1,counter)=diag(A,-1);
m=m-1;
A=A(1:m,1:m);
end
end
