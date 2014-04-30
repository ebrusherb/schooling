function h = H2norm(M,scaled)
L=lap(M);
s=size(L);
n=s(1);
if nargin==1
    scaled='additive';
end
if strcmp(scaled,'additive')==1
noise=diag(sum(transpose(M)));
else
noise=eye(n);
end
q=zeros(n-1,n);
q(1,1:2)=[1/sqrt(2) -1/sqrt(2)];
for j=2:(n-1)
   v=rand(1,n);
   v=v-sum(v)/n;
   for k=1:(j-1)
      v=v-(v*transpose(q(k,:)))*q(k,:); 
   end
   v=v/norm(v,2);
   q(j,:)=v;
end
Lbar=q*L*transpose(q);
noise=q*noise*transpose(q);
sigma = lyap(Lbar,-noise);
h = sqrt(trace(sigma));