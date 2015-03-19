strategy=13*ones(N,1);
positions=unifrnd(0,1,N,2);
d=squareform(pdist(positions));

M=zeros(N);
for ind=1:N
    [~, order]=sort(d(ind,:));
    neighbors=order(2:strategy(ind)+1);
    M(ind,neighbors)=1/strategy(ind);
end
M(1:N+1:end)=-1; %sets diagonal equal to -1
receivers=randsample(1,numsigs_permove,'true');
    
beta=zeros(N,1);
receiver=receivers(1);
allreceivers=d(receiver,:)<=radius;
beta(allreceivers)=b;
clc
t0=tic;q = expected_spin(M,T,beta);toc(t0)
t0=tic;q2 = expected_spin2(M,T,beta);toc(t0)
t0=tic;q_old = expected_spin_old(M,T,beta);toc(t0)
real(sigfig([q q2 q_old],10))


