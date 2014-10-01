function x = scores2(strategy,beta,probs)
x=cell(1,3);
M=makenet(strategy);
L=-lap(M);
s=size(M);
N=s(1);
dist=zeros(N);
noise=diag(sum(transpose(M)));
if nargin==2
    probs=1/N*ones(1,N); %probability of each node being the receiver
end
if var(sum(transpose(M)))==0 %if everyone's the same don't bother looping over all receivers
    L2=L;
    L2(1,1)=L2(1,1)+beta;
    sigma = lyap(L2,noise);
    for i=1:N
         dist(:,i) = circshift(reshape(diag(sigma),N,1),i-1);
    end
    dist=dist*diag(probs);
    x{1}=sum(diag(sigma))/N;
    x{2}=sum(sum(sigma))/N^2;
    x{3}=sum(dist,2);
end
if var(sum(transpose(M)))~=0
    consensus=zeros(1,N);
    suscept=zeros(1,N);
    for i=1:N
        L2=L;
        L2(i,i)=L2(i,i)+beta;
        sigma = lyap(L2,noise);
        dist(:,i) = diag(sigma);
        consensus(i)=sum(diag(sigma))/N;
        suscept(i)=sum(sum(sigma))/N^2;
    end
    dist=dist*diag(probs);
    consensus=consensus.*probs;
    suscept=suscept.*probs;
    x{1}=sum(consensus);
    x{2}=sum(suscept);
    x{3}=sum(dist,2);
end




