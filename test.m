close all
N=20;
strategy=2*ones(N,1);
b=1;
T=1;
radius=.1;

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

% gplotdc(M,positions)
% rainbow=colormap(jet(max(Ctotal)));
% hold on
% for i=1:max(Ctotal)
%     plot(positions(Ctotal==i,1),positions(Ctotal==i,2),'o','LineWidth',5,'Color',rainbow(i,:))
% end
% plot(positions(receiver,1),positions(receiver,2),'ok','MarkerSize',10,'LineWidth',3)
% colorbar

groups=netgroups(M);
ng=size(groups,1);

inliers=~logical((Ctotal==max(Ctotal)+1));
% inliers=~logical((Ctotal==5)+(Ctotal==2));
Mreduced=M(inliers,inliers);
betareduced=beta(inliers);
Nreduced=size(Mreduced,2);
positionsreduced=positions(inliers,:);

q=zeros(Nreduced-1,Nreduced);
q(1,1:2)=[1/sqrt(2) -1/sqrt(2)];
for j=2:(Nreduced-1)
   v=rand(1,Nreduced);
   v=v-sum(v)/Nreduced;
   for k=1:(j-1)
      v=v-(v*transpose(q(k,:)))*q(k,:); 
   end
   v=v/norm(v,2);
   q(j,:)=v;
end


A=Mreduced-diag(betareduced);
P=-q*A*transpose(q);
Pinv=inv(P);

[vecs,vals]=eig(P);
vals=diag(vals);
invvals=1./vals;

sigma1=q*transpose(A)*ones(Nreduced,1);
sigma2=q*betareduced;

zcorrs=zeros(Nreduced-1,Nreduced-1);
for i=1:(Nreduced-1)
    for j=i:(Nreduced-1)
        piece1=sum(vecs(i,:).*invvals'.*vecs(j,:));
        piece2=sum(sum((sigma1'*Pinv(:,i))*(sigma1'*Pinv(:,j))'))+sum(sum((sigma1'*Pinv(:,j))*(sigma1'*Pinv(:,i))'));
        piece3=sum(sum((sigma2'*Pinv(:,i))*(sigma2'*Pinv(:,j))'))+sum(sum((sigma2'*Pinv(:,j))*(sigma2'*Pinv(:,i))'));
        zcorrs(i,j)=piece1+piece2+piece3;
        zcorrs(j,i)=piece1+piece2+piece3;
    end
end

% corrs=transpose(q)*zcorrs*q;

cov=transpose(q)*zcorrs*q;
todivide=repmat(diag(cov),1,Nreduced).*repmat(diag(cov)',Nreduced,1);
todivide=power(todivide,.5);
corrs=cov./todivide;

Mbin=Mreduced;
Mbin(Mreduced==-1)=0;
Mbin(~~Mbin)=1;
g=sparse(Mbin);
[~,C]=graphconncomp(g,'Directed','true');
[~,Cun]=graphconncomp(g,'Directed','false');
[~,group]=sort(C);
[~,groupun]=sort(Cun);
[S,~]=graphconncomp(g);
[veca,vala]=eig(A);
vala=diag(vala);
[vecm,valm]=eig(Mreduced);
valm=diag(valm);
s=find(sigfig(valm,14)==0);
% real([sort(-vala) [NaN;sort(vals)] veca(group,s) C(group)'])
figure
gplotdc(Mreduced,positionsreduced)
rainbow=colormap(jet(max(C)));
hold on
for i=1:max(C)
    plot(positionsreduced(C==i,1),positionsreduced(C==i,2),'o','LineWidth',5,'Color',rainbow(i,:))
end
plot(positionsreduced(receiver,1),positionsreduced(receiver,2),'ok','MarkerSize',10,'LineWidth',3)
colorbar
set(gcf,'Position',[0 378 560 420])

figure
gplotdc(Mreduced,positionsreduced)
rainbow=colormap(jet(max(Cun)));
hold on
for i=1:max(Cun)
    plot(positionsreduced(Cun==i,1),positionsreduced(Cun==i,2),'o','LineWidth',5,'Color',rainbow(i,:))
end
plot(positionsreduced(receiver,1),positionsreduced(receiver,2),'ok','MarkerSize',10,'LineWidth',3)
colorbar
set(gcf,'Position',[750 378 560 420])

% [veca(group,s) C(group)' Cun(group)']
% [vecm(group,s) C(group)' Cun(group)']

if max(max(abs(imag(corrs))))==0
    figure
    imagesc(corrs(groupun,groupun))
    colorbar
    caxis manual
    caxis([-1 1])
    set(gcf,'Position',[750 78 560 420])
end
    
diag(corrs)
ng

%%
% inliers=~logical((Ctotal==max(Ctotal)+1));
% inliers=logical((Ctotal==2));
inliers=logical(~inliers);
Mreduced=M(inliers,inliers);
betareduced=beta(inliers);
Nreduced=size(Mreduced,2);
positionsreduced=positions(inliers,:);

q=zeros(Nreduced-1,Nreduced);
q(1,1:2)=[1/sqrt(2) -1/sqrt(2)];
for j=2:(Nreduced-1)
   v=rand(1,Nreduced);
   v=v-sum(v)/Nreduced;
   for k=1:(j-1)
      v=v-(v*transpose(q(k,:)))*q(k,:); 
   end
   v=v/norm(v,2);
   q(j,:)=v;
end


A=Mreduced-diag(betareduced);
P=-q*A*transpose(q);
Pinv=inv(P);

[vecs,vals]=eig(P);
vals=diag(vals);
invvals=1./vals;

sigma1=q*transpose(A)*ones(Nreduced,1);
sigma2=q*betareduced;

zcorrs=zeros(Nreduced-1,Nreduced-1);
for i=1:(Nreduced-1)
    for j=i:(Nreduced-1)
        piece1=sum(vecs(i,:).*invvals'.*vecs(j,:));
        piece2=sum(sum((sigma1'*Pinv(:,i))*(sigma1'*Pinv(:,j))'))+sum(sum((sigma1'*Pinv(:,j))*(sigma1'*Pinv(:,i))'));
        piece3=sum(sum((sigma2'*Pinv(:,i))*(sigma2'*Pinv(:,j))'))+sum(sum((sigma2'*Pinv(:,j))*(sigma2'*Pinv(:,i))'));
        zcorrs(i,j)=piece1+piece2+piece3;
        zcorrs(j,i)=piece1+piece2+piece3;
    end
end

% corrs=transpose(q)*zcorrs*q;

cov=transpose(q)*zcorrs*q;
todivide=repmat(diag(cov),1,Nreduced).*repmat(diag(cov)',Nreduced,1);
todivide=power(todivide,.5);
corrs=cov./todivide;

Mbin=Mreduced;
Mbin(Mreduced==-1)=0;
Mbin(~~Mbin)=1;
g=sparse(Mbin);
[~,C]=graphconncomp(g,'Directed','true');
[~,Cun]=graphconncomp(g,'Directed','false');
[~,group]=sort(C);
[~,groupun]=sort(Cun);
[S,~]=graphconncomp(g);
[veca,vala]=eig(A);
vala=diag(vala);
[vecm,valm]=eig(Mreduced);
valm=diag(valm);
s=find(sigfig(vala,14)==0);
% real([sort(-vala) [NaN;sort(vals)] veca(group,s) C(group)'])
figure
gplotdc(Mreduced,positionsreduced)
rainbow=colormap(jet(max(C)));
hold on
for i=1:max(C)
    plot(positionsreduced(C==i,1),positionsreduced(C==i,2),'o','LineWidth',5,'Color',rainbow(i,:))
end
plot(positionsreduced(receiver,1),positionsreduced(receiver,2),'ok','MarkerSize',10,'LineWidth',3)
colorbar
set(gcf,'Position',[0 378 560 420])

figure
gplotdc(Mreduced,positionsreduced)
rainbow=colormap(jet(max(Cun)));
hold on
for i=1:max(Cun)
    plot(positionsreduced(Cun==i,1),positionsreduced(Cun==i,2),'o','LineWidth',5,'Color',rainbow(i,:))
end
plot(positionsreduced(receiver,1),positionsreduced(receiver,2),'ok','MarkerSize',10,'LineWidth',3)
colorbar
set(gcf,'Position',[750 378 560 420])

% [veca(group,s) C(group)' Cun(group)']
[vecm(group,s) C(group)' Cun(group)']

if max(max(abs(imag(corrs))))==0
    figure
    imagesc(corrs(group,group))
    colorbar
    caxis manual
    caxis([-1 1])
    set(gcf,'Position',[750 78 560 420])
end
    