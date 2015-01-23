figure
% plot(eigcen,scores,'ob')
hold on


numrec=sum(keepreceivers,1);

mycolormap=cbrewer('qual','Set1',max(numrec));
for i=1:max(numrec)
    plot(eigcen(:,numrec==i),scores(:,numrec==i),'o','Color',mycolormap(i,:))
end
set(gca,'xlim',[.2 .3])
set(gca,'ylim',[.48 .6])
plot(eigcen(invader,:),scores(invader,:),'ok')

%%
n=5;
f=find(numrec==n);
plot(eigcen(:,f),scores(:,f),'o')
set(gca,'xlim',[.2 .3])
set(gca,'ylim',[.48 .6])

%%
figure
plot(repmat(eigval(:,f),N,1),eigcen(:,f),'o')
%%
f2=find(minscorer(invader,:)==.001);
f3=setdiff(1:1000,f2);
figure
% plot(eigcen,scores,'ob')
hold on

numrec=sum(keepreceivers,1);

mycolormap=cbrewer('qual','Set1',max(numrec));
for i=1:max(numrec)
    fnow=find(numrec==i);
    plot(eigcen(:,intersect(fnow,f2)),scores(:,intersect(fnow,f2)),'o','Color',mycolormap(i,:))
end
set(gca,'xlim',[.2 .3])
set(gca,'ylim',[.48 .6])
plot(eigcen(invader,f2),scores(invader,f2),'ok')

%%
k=k+1;
% plot(eigcen(:,f3(k)),scores(:,f3(k)),'o')
plot(eigcen(:,k),scores(:,k),'ok','MarkerSize',20)
hold on
plot(eigcen(invader,k),scores(invader,k),'or')
% hold off

%%
normsumtosignal=sumtosignal.*repmat(strategy',1,1000);
figure

hold on

mycolormap=cbrewer('qual','Set1',max(max(normsumtosignal)));
p=zeros(1,max(max(normsumtosignal)));
leglab=cell(1,max(max(normsumtosignal)));
% for i=1:max(max(normsumtosignal))
for i=1:3
    fnow=find(col(normsumtosignal)'==i);
    if ~isempty(fnow)
        p(i)=plot(eigcen(fnow),scores(fnow),'o','Color',mycolormap(i,:),'LineWidth',5);
        leglab{i}=num2str(i);
    end
end
legend(p(p>0),leglab(p>0))
% set(gca,'xlim',[.2 .3])
% set(gca,'ylim',[.48 .6])
plot(eigcen(invader,:),scores(invader,:),'ok')

%%
numrec=sum(keepreceivers,1);
normsumtosignal=sumtosignal.*repmat(strategy',1,1000);
resscores=scores([1:(invader-1) (invader+1):N],:);
resnormsumtosignal=normsumtosignal([1:(invader-1) (invader+1):N],:);
figure
hold on
mycolormap=cbrewer('qual','Set1',max(numrec));

for i=1:max(max(normsumtosignal))
    fnow=find(col(resnormsumtosignal)'==i);
    if ~isempty(fnow)
        if i<strat
%             plot(i,resscores(fnow),'o','Color',mycolormap(i,:),'LineWidth',5);
            plot(i,resscores(fnow),'or','LineWidth',3);
        else
%             plot(i-strat+1,resscores(fnow),'o','Color',mycolormap(i-strat+1,:),'LineWidth',5);
            plot(i-strat+1,resscores(fnow),'or','LineWidth',3);
        end
    end
end

for i=1:max(max(normsumtosignal))
    fnow=find(normsumtosignal(invader,:)==i);
    if ~isempty(fnow)
        if i<invaderstrat
            plot(i,scores(invader,fnow),'ok','LineWidth',1,'MarkerSize',10);
        else
            plot(i-invaderstrat+1,scores(invader,fnow),'ok','LineWidth',1,'MarkerSize',10);
        end
    end
end
