N=20;
radius=radvals(3);
T=Tvals(3);
nummoves=1000;
numsigs_permove=1;
b=1;
perm=randsample(20,20,'false');

strats=[16 19];
invaderstrats=[17 16];

maxrec=9;
invader_scores=cell(maxrec+1,2);
invader_counts=zeros(maxrec+1,2);
invader_minscorer=cell(maxrec+1,2);
resident_scores=cell(maxrec+1,2);
resident_counts=zeros(maxrec+1,2);
resident_minscorer=cell(maxrec+1,2);
ninv=zeros(3,2);
nres=zeros(3,2);

for k=1:2
    strat=strats(k);
    invader=10;
    invaderstrat=invaderstrats(k);
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

    numrec=sum(keepreceivers,1);
    normsumtosignal=sigfig(sumtosignal.*repmat(strategy',1,1000));
    normsumtosignal_res=normsumtosignal([1:(invader-1) (invader+1):end],:);
    scores_res=scores([1:(invader-1) (invader+1):end],:);
    minscorer_res=minscorer([1:(invader-1) (invader+1):end],:);
        
    for i=1:(maxrec+1)
        finv=find(normsumtosignal(invader,:)==i-1);
        invader_scores{i,k}=col(scores(invader,finv));
        invader_counts(i,k)=length(finv);
        invader_minscores{i,k}=col(minscorer(invader,finv));
        fres=find(normsumtosignal_res==i-1);
        resident_scores{i,k}=col(scores_res(fres));
        resident_counts(i,k)=length(fres);
        resident_minscores{i,k}=col(minscorer_res(fres));
    end
    
%     f=find(normsumtosignal(invader,:)==0);
%     v=minscorer(invader,f);
%     v=1./v/1000;
%     v=sigfig(v,2);
%     v(minscorer(invader,f)==0)=0;
%     [ninv(:,k),binv]=hist(v,0:12);
%     fres=find(normsumtosignal_res==0);
%     v=minscorer_res(fres);
%     v=1./v/1000;
%     v=sigfig(v,2);
%     v(minscorer_res(fres)==0)=0;
%     [nres(:,k),bres]=hist(v,0:12);
    
    vec=[0 .05 1];
    f=find(normsumtosignal(invader,:)==0);
    v=minscorer(invader,f)*1000;
    [ninv(:,k),binv]=hist(v,vec);
    fres=find(normsumtosignal_res==0);
    v=minscorer_res(fres)*1000;
    [nres(:,k),bres]=hist(v,vec);    

end

%%
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
subplot(2,2,1)
data=cell(maxrec+1,2);
for i=1:(maxrec+1)
    data{i,1}=resident_scores{i,1};
    data{i,2}=invader_scores{i,1};
    data{i,3}=resident_scores{i,2};
    data{i,4}=invader_scores{i,2};
end

xlab=cell(1,maxrec+1);
for i=1:(maxrec+1)
    xlab{i}=num2str(i-1);
end

cols=cbrewer('qual','Set1',4);
cols=horzcat(cols,200*ones(4,1)/255);
cols=cols';

Mlab={'Resident', 'Invader','Resident','Invader'};
p=multiple_boxplot(data,xlab,Mlab,cols,0);
xlabel('Number of informed neighbors')
ylabel('Score')

leg=legend(p,Mlab);
legend('boxoff')
set(leg,'Position',[.31 .65 .15 .07])

data2=cell(maxrec+1,2);
for i=1:(maxrec+1)
    data2{i,1}=resident_minscores{i,1}*1000;
    data2{i,2}=invader_minscores{i,1}*1000;
    data2{i,3}=resident_minscores{i,2}*1000;
    data2{i,4}=invader_minscores{i,2}*1000;
end

subplot(2,2,2)
Mlab={'Resident', 'Invader','Resident','Invader'};
p=multiple_boxplot(data2,xlab,Mlab,cols,0);
xlabel('Number of informed neighbors')
ylabel('Probability of being eaten')

% leg=legend(p,Mlab,'FontName','TimesNewRoman','FontSize',12);
% legend('boxoff')
% set(leg,'Position',[.3 .8 .15 .07])

subplot(2,2,3)
w=.3;
cols=cbrewer('qual','Set1',4);
    
bar((0:maxrec)-3/4*w,resident_counts(:,1)/sum(resident_counts(:,1)),w/2,'FaceColor',cols(1,:))
hold on
bar((0:maxrec)-w/4,invader_counts(:,1)/sum(invader_counts(:,1)),w/2,'FaceColor',cols(2,:));
bar((0:maxrec)+w/4,resident_counts(:,2)/sum(resident_counts(:,2)),w/2,'FaceColor',cols(3,:));
bar((0:maxrec)+3/4*w,invader_counts(:,2)/sum(invader_counts(:,2)),w/2,'FaceColor',cols(4,:));
% legend('Resident','Invader','Resident','Invader')
% legend('boxoff')
xlabel('Number of informed neighbors')
set(gca,'xtick',0:maxrec,'xticklabel',0:maxrec)

subplot(2,2,4)
w=.2;
bar([0 .5 1]-3/4*w,nres(:,1)/sum(nres(:,1)),w,'FaceColor',cols(1,:))
hold on
gb=bar([0 .5 1]-1/4*w,ninv(:,1)/sum(ninv(:,1)),w,'FaceColor',cols(2,:));
bar([0 .5 1]+1/4*w,nres(:,2)/sum(nres(:,2)),w,'FaceColor',cols(3,:))
bar([0 .5 1]+3/4*w,ninv(:,2)/sum(ninv(:,2)),w,'FaceColor',cols(4,:))
hold off
% legend('Resident','Invader','Resident','Invader')
% legend('boxoff')
xlabel('Probability of being eaten during signaling event')
set(gca,'xtick',[0 .5 1],'xticklabel',{'0','intermediate','1'})

% set(gcf,'PaperSize',[w h]);
% set(gcf,'PaperPosition',[-w*.55 0 w*2 h]);

filename=strcat('/Users/eleanorbrush/Desktop/','ESS_figures_test','.pdf');
% print(filename,'-dpdf','-r300');
%%
f=find(normsumtosignal(invader,:)==0);
v=minscorer(invader,f)*1000;
% vec=sort(1./(1:12)/1000-.5*(1/12/1000));
% vec=(0:.01:1)-.005;
vec=[0 .05 1];
[n,bins]=hist(v,vec);
bar([0 .5 1],n/sum(n));
set(gca,'xtick',[0 .5 1],'xticklabel',{'0','intermediate','1'})

%%
f1=find(minscorer(invader,:)==.001);
f2=find(normsumtosignal(invader,:)==0);
f3=intersect(f1,f2);

%%
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=4;
h=w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
xlab=cell(1,maxrec+1);
for i=1:(maxrec+1)
    xlab{i}=num2str(i-1);
end

cols=cbrewer('qual','Set1',4);
cols=horzcat(cols,200*ones(4,1)/255);
cols=cols';
Mlab={'Resident', 'Invader','Resident','Invader'};
p=multiple_boxplot(data2,xlab,Mlab,cols,0);
xlabel('Number of informed neighbors')
ylabel('Probability of being eaten')

leg=legend(p,Mlab,'FontName','TimesNewRoman','FontSize',12);
legend('boxoff')
set(leg,'Position',[.3 .8 .15 .07])

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[-w*.55 0 w*2 h]);

filename=strcat('/Users/eleanorbrush/Desktop/','prob_eaten_test','.pdf');
print(filename,'-dpdf','-r300');

%%
lu=13;
uvec=0:(lu-1);
totalprob=zeros(10,4);

for i=1:10
    for j=1:4
        for k=1:lu
        totalprob(i,j)=totalprob(i,j)+sum(data2{i,j}==uvec(k))*uvec(k)/1000;
        end
    end
end

cols=cbrewer('qual','Set1',4);
cols=horzcat(cols,200*ones(4,1)/255);
cols=cols';

hold on

p=zeros(1,4);
for j=1:4
    p(j)=plot(totalprob(:,j),'-o','Color',cols(1:3,j),'LineWidth',3);
end
Mlab={'Resident', 'Invader','Resident','Invader'};
leg=legend(p,Mlab);
legend('boxoff')
xlabel('Number of informed neighbors')
ylabel('Probability of being eaten')