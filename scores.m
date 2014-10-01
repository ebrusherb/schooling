function x = scores(strategy)
x=cell(1,3);
M=makenet(strategy);
L=lap(M);
s=size(L);
N=s(1);
noise=diag(sum(transpose(M)));
q=zeros(N-1,N);
q(1,1:2)=[1/sqrt(2) -1/sqrt(2)];
for j=2:(N-1)
   v=rand(1,N);
   v=v-sum(v)/N;
   for k=1:(j-1)
      v=v-(v*transpose(q(k,:)))*q(k,:); 
   end
   v=v/norm(v,2);
   q(j,:)=v;
end
Lbar=q*L*transpose(q);
noise=q*noise*transpose(q);
sigma = lyap(Lbar,-noise);
x{1}=trace(sigma)/N;
% x{2}=sum(sum(sigma));
dist=zeros(N);
range=1:N;
noise=diag(sum(transpose(M)));
probs=1/N*ones(1,N); %probability of each node being the receiver
if var(sum(transpose(M)))==0 %if everyone's the same don't bother looping over all receivers
    Li=L(range~=1,range~=1);
    sigma = lyap(Li,-noise(range~=1,range~=1));
    for i=1:N
         dist(range~=i,i) = circshift(reshape(diag(sigma),N-1,1),i-1);
    end
    dist=dist*diag(probs);
    x{2}=sum(sum(sigma))/N^2;
    x{3}=sum(dist,2);
end
if var(sum(transpose(M)))~=0
    suscept=zeros(1,N);
    for i=1:N
        Li=L(range~=i,range~=i);
        sigma = lyap(Li,-noise(range~=i,range~=i));
        dist(range~=i,i) = diag(sigma);
        suscept(i)=sum(sum(sigma))/N^2;
    end
    dist=dist*diag(probs);
    suscept=suscept.*probs;
    x{2}=sum(suscept);
    x{3}=sum(dist,2);
end



