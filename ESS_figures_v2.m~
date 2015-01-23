N=20;
radius=.1;
T=Tvals(3);
nummoves=5000;
numsigs_permove=1;
b=1;
perm=randsample(20,20,'false');

strats=[16 17];
invaderstrats=[17 16];

maxrec=zeros(1,2);
invader_scores=cell(N+1,2);
invader_counts=zeros(N+1,2);
invader_minscores=cell(N+1,2);
resident_scores=cell(N+1,2);
resident_counts=zeros(N+1,2);
resident_minscores=cell(N+1,2);
ninv=zeros(3,2);
nres=zeros(3,2);
probeaten=zeros(N,2);

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
    keepM=reshape(keepM,20,20,nummoves);
    keepreceivers=reshape(keepreceivers,[],nummoves);
    maxrec(k)=max(sum(keepreceivers,1));
    eigcen=reshape(eigcen,[],nummoves);
    eigval=reshape(eigval,[],nummoves);
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
    probeaten(:,k)=sum(minscorer,2);     

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
    normsumtosignal=sigfig(sumtosignal.*repmat(strategy',1,nummoves));
    normsumtosignal_res=normsumtosignal([1:(invader-1) (invader+1):end],:);
    scores_res=scores([1:(invader-1) (invader+1):end],:);
    minscorer_res=minscorer([1:(invader-1) (invader+1):end],:);
        
    for i=1:(N+1)
        finv=find(normsumtosignal(invader,:)==i-1);
        invader_scores{i,k}=col(scores(invader,finv));
        invader_counts(i,k)=length(finv);
        invader_minscores{i,k}=col(minscorer(invader,finv));
        fres=find(normsumtosignal_res==i-1);
        resident_scores{i,k}=col(scores_res(fres));
        resident_counts(i,k)=length(fres);
        resident_minscores{i,k}=col(minscorer_res(fres));
    end 

end


totalprob_invader=zeros(N+1,2);
totalprob_resident=zeros(N+1,2);

for i=1:(N+1)
    for k=1:2
        if invader_counts(i,k)>0
            totalprob_invader(i,k)=sum(invader_minscores{i,k})/.001/invader_counts(i,k);
        end
    end
end
for i=1:(N+1)
    for k=1:2
        if resident_counts(i,k)>0
            totalprob_resident(i,k)=sum(resident_minscores{i,k})/.001/resident_counts(i,k);
        end
    end
end


figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=.75*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);

cols=cbrewer('qual','Set1',4);
cols=horzcat(cols,200*ones(4,1)/255);
cols=cols';

for j=1:2
    subplot(1,2,j)
    hold on
    p=zeros(1,2);
    if j==1
        p(1)=plot(0:N,totalprob_resident(:,j),'-o','Color',cols(1:3,1),'LineWidth',3);
        p(2)=plot(0:N,totalprob_invader(:,j),'-o','Color',cols(1:3,2),'LineWidth',3);
    else
        p(1)=plot(0:N,totalprob_resident(:,j),'-o','Color',cols(1:3,2),'LineWidth',3);
        p(2)=plot(0:N,totalprob_invader(:,j),'-o','Color',cols(1:3,1),'LineWidth',3);
    end
    set(gca,'xlim',[0 max(maxrec)])
    xlabel('Number of informed neighbors')
    ylabel('Probability of being eaten')
    Mlab={strcat('Resident = ',num2str(strats(j))),strcat( 'Invader = ',num2str(invaderstrats(j)))};
    leg=legend(p,Mlab);
    legend('boxoff')
    v=get(leg,'Position');
    hold off
end

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','mutual_uninvasibility','.pdf');
print(filename,'-dpdf','-r300');
%%
figure
w=.3;
cols=cbrewer('qual','Set1',4);
    
bar((0:N)-3/4*w,resident_counts(:,1)/sum(resident_counts(:,1)),w/2,'FaceColor',cols(1,:))
hold on
bar((0:N)-w/4,invader_counts(:,1)/sum(invader_counts(:,1)),w/2,'FaceColor',cols(2,:));
bar((0:N)+w/4,resident_counts(:,2)/sum(resident_counts(:,2)),w/2,'FaceColor',cols(3,:));
bar((0:N)+3/4*w,invader_counts(:,2)/sum(invader_counts(:,2)),w/2,'FaceColor',cols(4,:));
% legend('Resident','Invader','Resident','Invader')
% legend('boxoff')
xlabel('Number of informed neighbors')
set(gca,'xtick',0:max(maxrec),'xticklabel',0:max(maxrec),'xlim',[-.5 12.5])