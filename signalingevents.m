function [meanscore scorevar probeaten probgettoeat]=signalingevents(strategy,numsigs_permove,nummoves,radius,b,T)
% strategy=randi([1 N/2-1],1,N)*2;
numsigs_tot=numsigs_permove*nummoves;
N=max(size(strategy));

scores=zeros(N,numsigs_tot);
H2norms=zeros(1,numsigs_tot);
susceptvals=zeros(1,numsigs_tot);

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
    j=1;
    while j<=numsigs_permove
        beta=zeros(N,1);
        receiver=receivers(j);
        allreceivers=d(receiver,:)<=radius;
        beta(allreceivers)=b;
        v=real(expected_spin(M,T,beta));
        
        if sum(abs(v)>1000)>=1
            j=j;
        elseif sum(isnan(v))>=1
            j=j;
        else  
            scores(:,(i-1)*numsigs_permove+j)=v;
%             [h,s]=networkprops(M,'uniform',beta);
%             H2norms((i-1)*numsigs_permove+j)=h;
%             susceptvals((i-1)*numsigs_permove+j)=s;
            j=j+1;
        end
    end
end
weird=max(abs(scores),[],1)>1000;
scores(:,weird)=[];
meanscore=mean(scores,2);
scorevar=var(scores,[],2);

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

% meanH2=mean(H2norms);
% meansuscept=mean(susceptvals);
% subplot(3,1,1)
% plot(strategy,(meanscore)','o')
% subplot(3,1,2)
% plot(strategy,scorevar','o')
% subplot(3,1,3)
% plot(H2norms)
end