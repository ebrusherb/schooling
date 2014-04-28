function p=payoffs(strategy,numsigs_permove,nummoves,radius,b)
numsigs_tot=numsigs_permove*nummoves;
N=max(size(strategy));

scores=zeros(N,numsigs_tot);

for i=1:nummoves
    positions=unifrnd(0,1,N,2);
    d=squareform(pdist(positions));

    M=zeros(N);
    for ind=1:N
        [~, order]=sort(d(ind,:));
        neighbors=order(2:strategy(ind)+1);
        M(ind,neighbors)=1/strategy(ind);
    end
    M(1:N+1:end)=-1; %sets diagonal equal to -1
    receivers=randsample(N,numsigs_permove,'true');
    for j=1:numsigs_permove
        beta=zeros(N,1);
        receiver=receivers(j);
        allreceivers=d(receiver,:)<=radius;
        beta(allreceivers)=b;
        scores(:,(i-1)*numsigs_permove+j)=expected_spin(M,1,beta);
        if sum(abs(scores(:,(i-1)*numsigs_permove+j))>1000)>=1
            break
        end
    end
end
p=mean(scores,2);
end