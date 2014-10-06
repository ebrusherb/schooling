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
covtotal=covtotal*sign(covtotal(receiver,receiver));
figure 
imagesc(covtotal)
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

ingroup_paths=paths(:,receiver)~=Inf;
ingroup_eig=groupvec==groupvec(receiver);
ingroup_cov=index(diag(covtotal)>0);

figure
gplotdc(M,positions)
% rainbow=colormap(jet(max(ng)));
% hold on
% for i=1:max(ng)
%     plot(positions(groupvec==i,1),positions(groupvec==i,2),'o','LineWidth',5,'Color',rainbow(i,:))
% end
rainbow=colormap(jet(max(Ctotal)));
hold on
for i=1:max(Ctotal)
    plot(positions(Ctotal==i,1),positions(Ctotal==i,2),'o','LineWidth',5,'Color',rainbow(i,:))
end
good=index(ingroup_paths);
l=length(good);
for i=1:l
    plot(positions(good(i),1),positions(good(i),2),'ok','LineWidth',3,'MarkerSize',10)
end
good=index(ingroup_cov);
l=length(good);
for i=1:l
    plot(positions(good(i),1),positions(good(i),2),'or','LineWidth',3,'MarkerSize',15)
end
colorbar
set(gcf,'Position',[400 378 400 300])

ingroup=ingroup_paths;
distvec_paths=[];
corrvec_paths=[];
if sum(ingroup)>1
Mgroup=M(ingroup,ingroup);
betagroup=beta(ingroup);
bgroup=b;
dgroup=d(ingroup,ingroup);
receivergroup=find(index(ingroup)==receiver);
[distvec,corrvec,corrlength,crosscorrs] = correlationlength_mat_single_v3(Mgroup,dgroup,bgroup,radius,receivergroup);
distvec_paths=distvec;
corrvec_paths=corrvec;
end

ingroup=ingroup_eig;
distvec_eig=[];
corrvec_eig=[];
if length(ingroup)>1
Mgroup=M(ingroup,ingroup);
betagroup=beta(ingroup);
bgroup=b;
dgroup=d(ingroup,ingroup);
receivergroup=find(index(ingroup)==receiver);
[distvec,corrvec,corrlength,crosscorrs] = correlationlength_mat_single_v3(Mgroup,dgroup,bgroup,radius,receivergroup);
distvec_eig=distvec;
corrvec_eig=corrvec;
end

ingroup=ingroup_cov;
distvec_cov=[];
corrvec_cov=[];
if length(ingroup)>1
Mgroup=M(ingroup,ingroup);
betagroup=beta(ingroup);
bgroup=b;
dgroup=d(ingroup,ingroup);
receivergroup=find(index(ingroup)==receiver);
[distvec,corrvec,corrlength,crosscorrs] = correlationlength_mat_single_v3(Mgroup,dgroup,bgroup,radius,receivergroup);
distvec_cov=distvec;
corrvec_cov=corrvec;
end

ingroup=index;
distvec_total=[];
corrvec_total=[];
if abs(sum(sign(diag(covtotal))))==N
Mgroup=M(ingroup,ingroup);
betagroup=beta(ingroup);
bgroup=b;
dgroup=d(ingroup,ingroup);
receivergroup=find(index(ingroup)==receiver);
[distvec,corrvec,corrlength,crosscorrs] = correlationlength_mat_single_v3(Mgroup,dgroup,bgroup,radius,receivergroup);
distvec_total=distvec;
corrvec_total=corrvec;
end

figure
hold on
plot(distvec_paths,corrvec_paths,'-ok','MarkerSize',10)
plot(distvec_eig,corrvec_eig,'-o')
plot(distvec_cov,corrvec_cov,'-or')
plot(distvec_total,corrvec_total,'-og')
legend('paths','eig','cov','total')
set(gcf,'Position',[800 378 400 300])