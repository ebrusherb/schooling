N=20;
strategy=2*ones(N,1);

deltad=.05;

dbins=(0:deltad:1.4)+deltad/2;
Nbins=length(dbins);

radvals=0:.1:1;
Nr=length(radvals);
corrlength_vec=zeros(1,Nr);
cov_mat=zeros(N,N,Nr);
corr_mat=zeros(N,N,Nr);
avgcorr_mat=zeros(Nbins,Nr);

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
if length(w)>1
    'prob'
end
Wtilde=W(:,setdiff(1:N,w));
Lambdatilde=Lambda(setdiff(1:N,w),setdiff(1:N,w));
Cov=Wtilde*inv(Lambdatilde)*transpose(Wtilde); %#ok<*MINV>
plot(d,-Cov,'o');set(gca,'ylim',[-.5 5])      
%%
%     noise=eye(N);
noise=diag(strategy);
    
receiver=1;
for j=1:Nr
    radius=radvals(j);
    allreceivers=d(receiver,:)<=radius;
    beta=zeros(N,1);
    beta(allreceivers)=b;
    B=diag(beta);
    Pf=P-B; 

    toinvert=-Pf-1/(b*sum(allreceivers))*B*ones(N,1)*ones(1,N)*B;
    [W,Lambda]=eig(toinvert);
    w=find(sigfig(diag(Lambda),13)==0);
    Cov_forced=W(:,[1:(w-1) (w+1):end])*inv(Lambda([1:(w-1) (w+1):end],[1:(w-1) (w+1):end]))*transpose(W(:,[1:(w-1) (w+1):end])); %#ok<MINV>
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
        w=intersect(w1,w2);
        avgcorrs_forced(i)=mean(Corr_vec(w));
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
        