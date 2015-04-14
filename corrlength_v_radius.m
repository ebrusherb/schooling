clc
close all
N=20;
strategy=13*ones(N,1);

deltad=.05;

dbins=(0:deltad:1.4)+deltad/2;
Nbins=length(dbins);
%%
radius=.1;
radvals=0:.1:1;
Nr=length(radvals);
corrlength_vec=zeros(1,Ns);
cov_mat=zeros(N,N,Ns);
corr_mat=zeros(N,N,Ns);
avgcorr_mat=zeros(Nbins,Ns);
numreceivers=zeros(1,Ns);

positions=unifrnd(0,1,N,2);
d=squareform(pdist(positions));
dvec_forced=col(d);
   
for j=1:Ns
    strategy=stratvals(j)*ones(N,1);
    A=zeros(N);
    Abar=zeros(N);
    for ind=1:N
        [~, order]=sort(d(ind,:));
        neighbors=order(2:strategy(ind)+1);
        A(ind,neighbors)=1;
        Abar(ind,neighbors)=1/strategy(ind);
    end

    S=.5*(Abar+transpose(Abar));
    P=S;
    P(1:N+1:end)=-sum(S,2);

    receiver=1;
    radius=radvals(j);
    allreceivers=d(receiver,:)<=radius;
    numreceivers(j+1)=sum(allreceivers);
    beta=zeros(N,1);
    beta(allreceivers)=b;
    B=diag(beta);
    Pf=P-B; 

    toinvert=-Pf-1/(b*sum(allreceivers))*B*ones(N,1)*ones(1,N)*B;
    [W,Lambda]=eig(toinvert);
    w=find(sigfig(diag(Lambda),13)==0);
    Cov_forced=W(:,setdiff(1:N,w))*inv(Lambda(setdiff(1:N,w),setdiff(1:N,w)))*transpose(W(:,setdiff(1:N,w))); 
    todivide=repmat(diag(Cov_forced),1,N).*repmat(reshape(diag(Cov_forced),1,[]),N,1);
    todivide=power(todivide,.5);
    Corr_forced=Cov_forced./todivide;
    Corr_vec=col(Corr_forced);
    cov_mat(:,:,j)=Cov_forced;
    corr_mat(:,:,j)=Corr_forced;

    avgcorrs_forced=zeros(1,Nbins);

    for i=1:Nbins
        w1=find(dvec_forced>=dbins(i)-deltad/2);
        w2=find(dvec_forced<dbins(i)+deltad/2);
        w3=intersect(w1,w2);
        avgcorrs_forced(i)=mean(Corr_vec(w3));
    end
    avgcorr_mat(:,j)=avgcorrs_forced;

    f=find(col(avgcorrs_forced)<0,1,'first');
    if isempty(f)
        corrlength_forced=max(dvec_forced);
    else
        corrlength_forced=interp1(avgcorrs_forced([f-1 f]),dbins([f-1 f]),0);
    end
    
    corrlength_vec(j)=corrlength_forced;
end

[s,o]=sort(d(receiver,:));
m=min(min(min(cov_mat)));
M=max(max(max(cov_mat)));
M=max(abs(m),M);

m2=min(min(min(corr_mat)));
M2=max(max(max(corr_mat)));
M2=max(abs(m2),M2);

figure
hold on
plot(dbins,avgcorr_mat(:,1),'Color',seqcols3(1,:))
for i=1:Nr
    plot(dbins,avgcorr_mat(:,i+1),'Color',seqcols3(i+1,:),'LineWidth',3)
end

legendlabs=cell(Nr+1,1);
for i=1:Nr
    legendlabs{i+1}=num2str(radvals(i));
end
legendlabs{1}='none';

legend(legendlabs)

for i=1:Nr
    plot(corrlength_vec(i),0,'o','Color',seqcols3(i,:),'LineWidth',5)
end

%%
stratvals=1:(N-1);
Ns=length(stratvals);
corrlength_vec=zeros(1,Nr+1);
cov_mat=zeros(N,N,Nr+1);
corr_mat=zeros(N,N,Nr+1);
avgcorr_mat=zeros(Nbins,Nr+1);
numreceivers=zeros(1,Nr+1);

positions=unifrnd(0,1,N,2);
d=squareform(pdist(positions));
dvec_forced=col(d);
    
A=zeros(N);
Abar=zeros(N);
for ind=1:N
    [~, order]=sort(d(ind,:));
    neighbors=order(2:strategy(ind)+1);
    A(ind,neighbors)=1;
    Abar(ind,neighbors)=1/strategy(ind);
end

S=.5*(Abar+transpose(Abar));
P=S;
P(1:N+1:end)=-sum(S,2);
[W,Lambda]=eig(P);
w=find(sigfig(diag(Lambda),13)==0);

Wtilde=W(:,setdiff(1:N,w));
Lambdatilde=Lambda(setdiff(1:N,w),setdiff(1:N,w));
Cov=Wtilde*inv(Lambdatilde)*transpose(Wtilde); %#ok<*MINV>
q=W(:,w)-1/(N)*repmat(sum(W(:,w),1),N,1).*ones(N,length(w));
[s,o]=sort(q(:,1));

todivide=repmat(diag(-Cov),1,N).*repmat(reshape(diag(-Cov),1,[]),N,1);
todivide=power(todivide,.5);
Corr=-Cov./todivide;
Corr_vec=col(Corr);
cov_mat(:,:,1)=-Cov;
corr_mat(:,:,1)=Corr;

avgcorrs_forced=zeros(1,Nbins);

for i=1:Nbins
    w1=find(dvec_forced>=dbins(i)-deltad/2);
    w2=find(dvec_forced<dbins(i)+deltad/2);
    w3=intersect(w1,w2);
    avgcorrs_forced(i)=mean(Corr_vec(w3));
end
avgcorr_mat(:,1)=avgcorrs_forced;

f=find(col(avgcorrs_forced)<0,1,'first');
if isempty(f)
    corrlength_forced=max(dvec_forced);
else
    corrlength_forced=interp1(avgcorrs_forced([f-1 f]),dbins([f-1 f]),0);
end

corrlength_vec(1)=corrlength_forced;

noise=diag(strategy);
    
receiver=1;
for j=1:Nr
    radius=radvals(j);
    allreceivers=d(receiver,:)<=radius;
    numreceivers(j+1)=sum(allreceivers);
    beta=zeros(N,1);
    beta(allreceivers)=b;
    B=diag(beta);
    Pf=P-B; 

    toinvert=-Pf-1/(b*sum(allreceivers))*B*ones(N,1)*ones(1,N)*B;
    [W,Lambda]=eig(toinvert);
    w=find(sigfig(diag(Lambda),13)==0);
    Cov_forced=W(:,setdiff(1:N,w))*inv(Lambda(setdiff(1:N,w),setdiff(1:N,w)))*transpose(W(:,setdiff(1:N,w))); %#ok<MINV>
    todivide=repmat(diag(Cov_forced),1,N).*repmat(reshape(diag(Cov_forced),1,[]),N,1);
    todivide=power(todivide,.5);
    Corr_forced=Cov_forced./todivide;
    Corr_vec=col(Corr_forced);
    cov_mat(:,:,j+1)=Cov_forced;
    corr_mat(:,:,j+1)=Corr_forced;

    avgcorrs_forced=zeros(1,Nbins);

    for i=1:Nbins
        w1=find(dvec_forced>=dbins(i)-deltad/2);
        w2=find(dvec_forced<dbins(i)+deltad/2);
        w3=intersect(w1,w2);
        avgcorrs_forced(i)=mean(Corr_vec(w3));
    end
    avgcorr_mat(:,j+1)=avgcorrs_forced;

    f=find(col(avgcorrs_forced)<0,1,'first');
    if isempty(f)
        corrlength_forced=max(dvec_forced);
    else
        corrlength_forced=interp1(avgcorrs_forced([f-1 f]),dbins([f-1 f]),0);
    end
    
    corrlength_vec(j+1)=corrlength_forced;
end

[s,o]=sort(d(receiver,:));
m=min(min(min(cov_mat)));
M=max(max(max(cov_mat)));
M=max(abs(m),M);

m2=min(min(min(corr_mat)));
M2=max(max(max(corr_mat)));
M2=max(abs(m2),M2);

figure
hold on
plot(dbins,avgcorr_mat(:,1),'Color',seqcols3(1,:))
for i=1:Nr
    plot(dbins,avgcorr_mat(:,i+1),'Color',seqcols3(i+1,:),'LineWidth',3)
end

legendlabs=cell(Nr+1,1);
for i=1:Nr
    legendlabs{i+1}=num2str(radvals(i));
end
legendlabs{1}='none';

legend(legendlabs)

for i=1:Nr
    plot(corrlength_vec(i),0,'o','Color',seqcols3(i,:),'LineWidth',5)
end
%%
figure
plot(corrlength_vec,'-o')

[s2,o2]=sort(d(1,:));
figure
imagesc(corr_mat(o2,o2,1))
colorbar
colormap(divcols)
caxis manual;
caxis([-M2 M2])

figure
k=1;
f2=find(d(1,o)>=corrlength_vec(k),1,'first');imagesc(corr_mat(o2,o2,k));colormap(divcols);colorbar;caxis manual;caxis([-M2 M2]);hold on;plot(k*ones(2,1),get(gca,'ylim'));plot(get(gca,'xlim'),k*ones(2,1));plot(f2*ones(2,1),get(gca,'ylim'),'r');plot(get(gca,'xlim'),f2*ones(2,1),'r');colorbar;caxis manual;caxis([-M2 M2])
        