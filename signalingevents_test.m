N=20;
radius=radvals(3);
T=Tvals(3);
nummoves=1000;
numsigs_permove=1;
b=1;
perm=randsample(20,20,'false');

strat=16;
invader=10;
invaderstrat=17;
strategy=[strat*ones(1,N)];
strategy(invader)=invaderstrat;

numsigs_tot=numsigs_permove*nummoves;
N=max(size(strategy));
indvec=1:N;

scores=zeros(N,numsigs_permove,nummoves);
keepreceivers=zeros(N,numsigs_permove,nummoves);
keepM=zeros(N,N,numsigs_permove,nummoves);
eigcen=zeros(N,numsigs_permove,nummoves);
eigval=zeros(numsigs_permove,nummoves);

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
        
        A=M-diag(beta);
        [vecs,vals]=eig(A);
        V=vecs;
        vals=diag(vals);
        nonzero=abs(vals)>0.00001;
        nonzerovals=vals(nonzero);
        l=1./nonzerovals.*(exp(nonzerovals*T)-1);
        Linhom=zeros(1,N);
        Linhom(nonzero)=l;
        Linhom(~nonzero)=T;
        v=V(:,Linhom==max(Linhom));
        eigcen(:,j,i)=v/sign(v(1));
        eigval(j,i)=max(Linhom);
    end
    
end
keepM=reshape(keepM,20,20,1000);
keepreceivers=reshape(keepreceivers,[],1000);
eigcen=reshape(eigcen,[],1000);
eigval=reshape(eigval,[],1000);
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
set(gcf,'Position',[10 300 700 500])
subplot(2,2,1)
stdev=std(scores,[],2);
minsc=min(scores,[],2);
maxsc=max(scores,[],2);
p1=plot(mean(scores,2),'green');
hold on
p2=plot(mean(scores,2)-stdev);
plot(mean(scores,2)+stdev)
p3=plot(minsc,'red');
plot(maxsc,'red')
hold off
legend([p1 p2 p3],{'Mean','stdev','min/max'})

subplot(2,2,2)
hist(rows,20)
xlabel('Individual')
ylabel('# of times ind is a min scorer')

% hist(minscorer,20)

% subplot(2,3,3)
% plot(mean(d))


payingattention=keepM;
payingattention(keepM==-1)=1;
tosignal=payingattention;
tonaive=payingattention;
for i=1:nummoves
    tosignal(:,:,i)=payingattention(:,:,i).*repmat(col(keepreceivers(:,i))',N,1);
    tonaive(:,:,i)=payingattention(:,:,i).*repmat(col(~keepreceivers(:,i))',N,1);
end
sumtosignal=reshape(sum(tosignal,2),N,nummoves);
sumtonaive=reshape(sum(tonaive,2),N,nummoves);

subplot(2,2,3)
hist(minscorer(invader,:)*1000,[-.05 :.1:1.1])
set(gca,'ylim',[0 1000])
xlabel('Prob eaten at each signaling event')
title('Invader')
subplot(2,2,4)
hist(minscorer(2,:)*1000,[-.05 :.1:1.1])
set(gca,'ylim',[0 1000])
xlabel('Prob eaten at each signaling event')
title('Resident')
% subplot(2,3,6)
% hist(rows,20)

% figure
% set(gcf,'Position',[800 300 700 500])
% subplot(2,2,1)
% hold on
% for i=1:N
%     p1=plot(strategy(i)*sumtosignal(i,:),minscorer(i,:),'o','LineWidth',5);
% end
% p2=plot(strategy(i)*sumtosignal(invader,:),minscorer(invader,:),'or','LineWidth',5,'MarkerSize',12);
% xlabel('Sum of weights to informed individuals')
% ylabel('Prob eaten at signaling event')
% legend([p1 p2],{'Resident','Invader'})
% subplot(2,2,2)
% hold on
% for i=1:N
%     plot(strategy(i)*sumtonaive(i,:),minscorer(i,:),'o','LineWidth',5)
% end
% plot(strategy(i)*sumtonaive(invader,:),minscorer(invader,:),'or','LineWidth',3,'MarkerSize',12)
% xlabel('Sum of weights to naive individuals')
% ylabel('Prob eaten at signaling event')
% subplot(2,2,3)
% [n1,x1]=hist(sumtosignal(invader,:),.1:.2:1.9);
% hold on
% bar(x1,n1,'r')
% [n2,x2]=hist(sumtonaive(invader,:),.1:.2:1.9);
% % bar(x2,n2,'b')
% set(gca,'xlim',[0 2],'ylim',[0 700])
% xlabel('Sum of weights to informed individuals')
% title('Invader')
% subplot(2,2,4)
% [n11,x11]=hist(sumtosignal(2,:),.1:.2:1.9);
% hold on
% bar(x11,n11,'r')
% [n21,x21]=hist(sumtonaive(2,:),.1:.2:1.9);
% % bar(x21,n21,'b')
% set(gca,'xlim',[0 2],'ylim',[0 700])
% xlabel('Sum of weights to informed individuals')
% title('Resident')

% figure
% hold on
% numrec=sum(keepreceivers,1);
% [~,numclasses]=sort(numrec);
% 
% mycolormap=cbrewer('qual','Set1',max(numrec));
% for i=1:max(numrec)
%     plot(eigcen(:,numrec==i),scores(:,numrec==i),'o','Color',mycolormap(i,:))
% end
% set(gca,'xlim',[.2 .3])
% set(gca,'ylim',[.48 .6])
% plot(eigcen(invader,:),scores(invader,:),'ok')

%%
normsumtosignal=sigfig(sumtosignal.*repmat(strategy',1,1000));
normsumtosignal_res=normsumtosignal([1:(invader-1) (invader+1):end],:);
scores_res=scores([1:(invader-1) (invader+1):end],:);
numrec=sum(keepreceivers,1);
maxrec=max(numrec);
invader_scores=cell(maxrec,1);
invader_counts=zeros(maxrec,1);
resident_scores=cell(maxrec,1);
resident_counts=zeros(maxrec,1);

for i=1:maxrec
    finv=find(normsumtosignal(invader,:)==i);
    invader_scores{i}=col(scores(invader,finv));
    invader_counts(i)=length(finv);
    fres=find(normsumtosignal_res==i);
    resident_scores{i}=col(scores_res(fres));
    resident_counts(i)=length(fres);
end

toboxplot_invader=zeros(1,sum(invader_counts));
labels_invader=zeros(1,sum(invader_counts));
toboxplot_resident=zeros(1,sum(resident_counts));
labels_resident=zeros(1,sum(resident_counts));

count=0;
for i=1:maxrec
    toboxplot_invader(count+(1:invader_counts(i)))=invader_scores{i};
    labels_invader(count+(1:invader_counts(i)))=i;
    count=count+invader_counts(i);
end

count=0;
for i=1:maxrec
    toboxplot_resident(count+(1:resident_counts(i)))=resident_scores{i};
    labels_resident(count+(1:resident_counts(i)))=i;
    count=count+resident_counts(i);
end

w=.3;

% figure
% gb=bar((1:maxrec)-w/2,invader_counts/sum(invader_counts),w);
% hold on
% bar((1:maxrec)+w/2,resident_counts/sum(resident_counts),w,'FaceColor','red')
% legend('Invader','Resident')
% xlabel('Number of informed neighbors')
% legend('boxoff')
% 
% figure
% hold on
% 
data=horzcat(resident_scores,invader_scores);

xlab=cell(1,maxrec);
for i=1:maxrec
    xlab{i}=num2str(i);
end

cols=[102,255,255, 200;
    51,153,255, 200];
cols=cols/255;
cols=cols';

Mlab={'Resident', 'Invader'};
p=multiple_boxplot(data,xlab,Mlab,cols,0);
xlabel('Number of informed neighbors')
ylabel('Score')

leg=legend(p,Mlab,'FontName','TimesNewRoman','FontSize',12);
legend('boxoff')
set(leg,'Position',[.3 .8 .15 .07])