%classical gram Schmidt code
function [q,r] = ndcgs(A)
[m,n]=size(A);
for j=1:n
v(:,j)=A(:,j);
for i=1:j-1
r(i,j)=transpose(q(:,i))*A(:,j);
v(:,j)=v(:,j)-r(i,j).*q(:,i);
end
r(j,j)=norm(v(:,j)); %2norm for a matrix is just norm;
q(:,j)=v(:,j)./r(j,j);
end
end
%modified gram Schmidt code
function [q,r] = ndmgs(A)
[m,n]=size(A);
for i=1:n
v(:,i)=A(:,i);
end
for i=1:n
r(i,i)=norm(v(:,i));
q(:,i)=v(:,i)/r(i,i);
for j=i+1:n
r(i,j)=transpose(q(:,i))*v(:,j);
v(:,j)=v(:,j)-r(i,j).*q(:,i);
end
end
end
