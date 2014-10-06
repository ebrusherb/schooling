close all
N=20;
strategy=2*ones(N,1);
b=1;
T=1;
radius=.1;
index=1:N;

positions=unifrnd(0,1,N,2);
d=squareform(pdist(positions));

M=zeros(N);
for ind=1:N
    [~, order]=sort(d(ind,:));
    neighbors=order(2:(strategy(ind)+1));
    M(ind,neighbors)=1/strategy(ind);
end
M(1:N+1:end)=-1;

beta=zeros(N,1);
receiver=1;
allreceivers=d(receiver,:)<=radius;
beta(allreceivers)=b;

Mbin=M;
Mbin(M==-1)=0;
Mbin(~~Mbin)=1;
g=sparse(Mbin);
[~,Ctotal]=graphconncomp(g,'Directed','false');
[paths] = graphallshortestpaths(g);

covtotal=paircovariance(M,beta);
figure 
imagesc(corrtotal)
colorbar
caxis manual
caxis([-1 1]);
set(gcf,'Position',[00 378 400 300])

groups=netgroups(M);
ng=size(groups,1);
groupvec=zeros(N,1);
for i=1:ng
    for j=1:length(groups{i})
        groupvec(Ctotal==groups{i}(j))=i;
    end
end

figure
gplotdc(M,positions)
rainbow=colormap(jet(max(ng)));
hold on
for i=1:max(ng)
    plot(positions(groupvec==i,1),positions(groupvec==i,2),'o','LineWidth',5,'Color',rainbow(i,:))
end
good=index(diag(covtotal)>0);
l=length(good);
for i=1:l
    plot(positions(good(i),1),positions(good(i),2),'ok','LineWidth',3,'MarkerSize',10)
end
% plot(positions(receiver,1),positions(receiver,2),'ok','MarkerSize',10,'LineWidth',3)
colorbar
set(gcf,'Position',[400 378 400 300])