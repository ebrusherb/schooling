function L=lap(M)
s=size(M);
n=s(1);
L=eye(n);
for i=1:n
    neighbors=M(i,:)>0;
    deg=sum(neighbors);
    L(i,neighbors)=L(i,neighbors)-1/deg;
end
