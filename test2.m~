close all
N=20;
strategy=3*ones(N,1);
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

groups=netgroups(M);
ng=size(groups,1);
groupvec=zeros(N,1);
for i=1:ng
    for j=1:length(groups{i})
        groupvec(Ctotal==groups{i}(j))=i;
    end
end

gplotdc(M,positions)
rainbow=colormap(jet(max(ng)));
hold on
for i=1:max(ng)
    plot(positions(groupvec==i,1),positions(groupvec==i,2),'o','LineWidth',5,'Color',rainbow(i,:))
end
plot(positions(receiver,1),positions(receiver,2),'ok','MarkerSize',10,'LineWidth',3)
colorbar

corrcell=cell(ng,1);
distveccell=cell(ng,1);
corrveccell=cell(ng,1);
corrlengthcell=cell(ng,1);
for i=1:ng
    Mgroup=M(groupvec==i,groupvec==i);
    betagroup=beta(groupvec==i);
    dgroup=d(groupvec==i,groupvec==i);
    bgroup=groupvec(receiver)==i;
    if groupvec(receiver)==i
        receivergroup=find(index(groupvec==i)==receiver);
    else receivergroup=1;
    end
    if size(Mgroup,1)>1
        corrcell{i}=paircorrelations(Mgroup,betagroup);
        [distvec,corrvec,corrlength,crosscorrs] = correlationlength_mat_single_v3(Mgroup,dgroup,bgroup,radius,receivergroup);
        distveccell{i}=distvec;
        corrveccell{i}=corrvec;
        corrlengthcell{i}=corrlength;
    end
end

groupnow=groupvec(receiver);
figure
plot(distveccell{groupnow},corrveccell{groupnow},'-o')
    