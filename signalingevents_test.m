N=20;
radius=1.5;
T=1;
nummoves=1000;
numsigs_permove=1;
b=1;
perm=randsample(20,20,'false');

strat=16;
invader=1;
invaderstrat=17;
strategy=[strat*ones(1,N)];
strategy(invader)=invaderstrat;
numsigs_tot=numsigs_permove*nummoves;
N=max(size(strategy));
indvec=1:N;

scores=zeros(N,numsigs_permove,nummoves);
keepreceivers=zeros(N,numsigs_permove,nummoves);
keepM=zeros(N,N,numsigs_permove,nummoves);

% parfor i=1:nummoves
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
%     M=M(perm,perm);
    receivers=randsample(N,numsigs_permove,'true');
    for j=1:numsigs_permove
        beta=zeros(N,1);
        receiver=receivers(j);
        allreceivers=d(receiver,:)<=radius;
        beta(allreceivers)=b;
        v=real(expected_spin(M,T,beta));
       
        scores(:,j,i)=v;
        keepreceivers(:,j,i)=(allreceivers);
        keepM(:,:,j,i)=M;
    end
end
keepM=reshape(keepM,20,20,1000);
keepreceivers=reshape(keepreceivers,[],1000);
scores=reshape(scores,N,[]);
% scores=scores(end:-1:1,:);
scores=sigfig(scores,4);
[c,o]=sort(scores,1);

% minscorer=o(1,:);


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

figure
subplot(2,2,1)
stdev=std(scores,[],2);
minsc=min(scores,[],2);
maxsc=max(scores,[],2);
plot(mean(scores,2),'green')
hold on
plot(mean(scores,2)-stdev)
plot(mean(scores,2)+stdev)
plot(minsc,'red')
plot(maxsc,'red')
hold off

subplot(2,2,2)
hist(rows,20)
% hist(minscorer,20)

% subplot(2,3,3)
% plot(mean(d))


payingattention=keepM;
payingattention(keepM==-1)=1;
tosignal=payingattention;
tonaive=payingattention;
for i=1:nummoves
    tosignal(:,:,i)=payingattention(:,:,i).*repmat(col(keepreceivers(:,i)),1,N);
    tonaive(:,:,i)=payingattention(:,:,i).*repmat(col(~keepreceivers(:,i)),1,N);
end
sumtosignal=reshape(sum(tosignal,2),N,nummoves);
sumtonaive=reshape(sum(tonaive,2),N,nummoves);

subplot(2,2,3)
hist(minscorer(invader,:)*1000,[-.05 :.1:1.1])
set(gca,'ylim',[0 1000])
subplot(2,2,4)
hist(minscorer(2,:)*1000,[-.05 :.1:1.1])
set(gca,'ylim',[0 1000])
% subplot(2,3,6)
% hist(rows,20)

figure
subplot(1,2,1)
hold on
for i=1:N
    plot(sumtosignal(i,:),minscorer(i,:),'o')
end
plot(sumtosignal(invader,:),minscorer(invader,:),'or')
subplot(1,2,2)
hold on
for i=1:N
    plot(sumtonaive(i,:),minscorer(i,:),'o')
end
plot(sumtonaive(invader,:),minscorer(invader,:),'or')


set(gcf,'Position',[50 50 600 1000])
% [probeaten, probgettoeat]=signalingevents(strategy,numsigs_permove,nummoves,radius,b,T);