function [probeaten, probgettoeat, meanlambda, meanH2, meancorrlength]=signalingevents_wgroupprops(strategy,numsigs_permove,nummoves,radius,b,T)
numsigs_tot=numsigs_permove*nummoves;
N=max(size(strategy));

scores=zeros(N,numsigs_permove,nummoves);
lambdavals=zeros(1,nummoves);
H2norms=zeros(1,nummoves);
corrlengths=zeros(numsigs_permove,nummoves);

parfor pari=1:nummoves
% for pari=1:nummoves
    
    positions=unifrnd(0,1,N,2);
    d=squareform(pdist(positions));

    M=zeros(N);
    for ind=1:N
        [~, order]=sort(d(ind,:));
        neighbors=order(2:(strategy(ind)+1));
        M(ind,neighbors)=1/strategy(ind);
    end
    M(1:N+1:end)=-1; %sets diagonal equal to -1
    
    L=lap(M);
    [~,vals]=eig(L);
    vals=diag(vals);
    vals=vals(2:length(vals));
    lambda=min(real(vals));
    lambdavals(pari)=lambda;
    h=H2norm(M,'additive');
    H2norms(pari)=h;   
    
    receivers=randsample(N,numsigs_permove,'true');
    for j=1:numsigs_permove
        beta=zeros(N,1);
        receiver=receivers(j);
        allreceivers=d(receiver,:)<=radius;
        beta(allreceivers)=b;
        v=real(expected_spin(M,T,beta));

        scores(:,j,pari)=v;

        [~,~,l,~]=correlationlength_mat_single_v3(M,d,b,radius,receiver);
        corrlengths(j,pari)=l;
    end
end

meanlambda=mean(lambdavals);
meanH2=mean(H2norms);
meancorrlength=mean(col(corrlengths));

scores=reshape(scores,N,[]);

[minvals,~]=min(scores);
minmat=repmat(minvals,N,1);
[rows,cols]=find(abs(scores-minmat)<0.00001);
minscorer=zeros(N,size(scores,2));

for i=1:size(scores,2)
    look=find(cols==i);
    minscorer(rows(look),i)=1/size(look,1)/numsigs_tot;
%     minscorer(rows(cols==i),i)=1/numsigs_tot;
end
probeaten=sum(minscorer,2);

[maxvals,~]=max(scores);
maxmat=repmat(maxvals,N,1);
[rows,cols]=find(abs(scores-maxmat)<0.00001);
maxscorer=zeros(N,size(scores,2));

for i=1:size(scores,2)
%     look=find(cols==i);
%     maxscorer(rows(look),i)=1/size(look,1)/numsigs_tot;
    maxscorer(rows(cols==i),i)=1/numsigs_tot;
end
probgettoeat=sum(maxscorer,2);
end