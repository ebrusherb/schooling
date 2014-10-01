function v=dynamics(x,M)
s=size(x);
n=s(1);
v=zeros(n,1);
for i=1:n
    neighbors=M(i,:)>0;
    deg=sum(neighbors);
   v(i,1)=x(i,1)+sum(x(neighbors,1)-x(i,1))/deg;
end
