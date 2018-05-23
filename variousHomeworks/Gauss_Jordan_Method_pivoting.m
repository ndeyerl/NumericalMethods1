function x=Gauss_Jordan_Method_pivoting(A,b) %Solve the Linear System Ax=b for x
if size(A,2)~=size(A,1);%if the matrix is not square
disp('Not a linear system of n equations with n unknowns ');   
elseif det(A)==0%if the matrix is not invertible
disp('A is not invertible ');
elseif size(b,1)~=size(A,2)%if the system is inconsistent
disp('dimensions issue'); 
else
n=size(A,2);%number of columns in A
c=zeros(n,1);% vector c initialized with zeros
x=zeros(n,1);% vector x initialized with zeros
[P,L,U]=PA_equal_LU(A);% PA=LU factorization of A
PB=P*b;
    for i=1:n%first backsubstitution, solve L*c=P*b for c 
       c(i)=PB(i)-L(i,1:i-1)*c(1:i-1); 
    end
    for i=n:-1:1%first backsubstitution, solve U*x=c for x
       x(i)=(c(i)-U(i,i+1:n)*x(i+1:n))/U(i,i); 
    end
end
end