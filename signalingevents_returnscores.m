function scores=signalingevents_returnscores(strategy,numsigs_permove,nummoves,radius,b,T)
N=max(size(strategy));

numsigs_tot=numsigs_permove*nummoves;


scores=zeros(N,numsigs_permove,nummoves);

% parfor i=1:nummoves
for i=1:nummoves
    positions=unifrnd(0,1,N,2);
    d=squareform(pdist(positions));

    M=zeros(N);
    for ind=1:N
        [~, order]=sort(d(ind,:));
        neighbors=order(2:strategy(ind)+1);
        M(ind,neighbors)=1/strategy(ind);
        M(ind,ind)=-sum(M(ind,neighbors));
    end
%     M(1:N+1:end)=-1; %sets diagonal equal to -1
%     for ind=1:N
%         if strategy(ind)==0
%             M(ind,ind)=0;
%         end
%     end
    receivers=randsample(N,numsigs_permove,'true');
    for j=1:numsigs_permove
        beta=zeros(N,1);
        receiver=receivers(j);
        allreceivers=d(receiver,:)<=radius;
        beta(allreceivers)=b;
        v=real(expected_spin(M,T,beta));
       
        scores(:,j,i)=v;
    end
end

scores=reshape(scores,N,[]);
