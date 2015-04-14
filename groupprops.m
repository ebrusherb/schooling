% function [meanH2, meanH2_forced, corrlength, corrlength_forced, disconnectedcount, meanrho, meanrho_forced]=groupprops(strategy,numsigs_permove,nummoves,radius,b,T) %#ok<INUSD>

% N=max(size(strategy));

dvec=[];
corrvec=[];
dvec_forced=[];
corrvec_forced=[];
disconnectedcount=0;
nooneislistening=0;
H2vec=zeros(nummoves,1);
H2vec_forced=zeros(nummoves,numsigs_permove);

Q=zeros(N-1,N);
Q(1,1:2)=[1/sqrt(2) -1/sqrt(2)];
for j=2:(N-1)
   v=rand(1,N);
   v=v-sum(v)/N;
   for k=1:(j-1)
      v=v-(v*transpose(Q(k,:)))*Q(k,:); 
   end
   v=v/norm(v,2);
   Q(j,:)=v;
end
Qfull=[Q; 1/sqrt(N)*ones(1,N)];


for i=1:nummoves   
    positions=unifrnd(0,1,N,2);
    d=squareform(pdist(positions));

    A=zeros(N);
    Abar=zeros(N);
    for ind=1:N
        [~, order]=sort(d(ind,:));
        neighbors=order(2:strategy(ind)+1);
        A(ind,neighbors)=1;
        Abar(ind,neighbors)=1/strategy(ind);
    end
  
    L=Abar-eye(N);
    Ltilde=Q*L*transpose(Q);
    [~,vals]=eig(Ltilde);
    problem=find(sigfig(diag(vals),13)==0,1);
    
%     noise=eye(N);
    strategynoise=strategy;
    strategynoise(strategy==0)=1;
    noise=diag(strategynoise);
    
    S=.5*(Abar+transpose(Abar));
    P=S;
    P(1:N+1:end)=-sum(S,2);
    [W,Lambda]=eig(P);
    w=find(sigfig(diag(Lambda),13)==0);
    
    Wtilde=W(:,setdiff(1:N,w));
    Lambdatilde=Lambda(setdiff(1:N,w),setdiff(1:N,w));
    Cov=Wtilde*inv(Lambdatilde)*transpose(Wtilde); %#ok<*MINV>

    todivide=repmat(diag(-Cov),1,N).*repmat(reshape(diag(-Cov),1,[]),N,1);
    todivide(sigfig(Cov,13)==0)=1;
    todivide=power(todivide,.5);
    Corr=-Cov./todivide;
    dvec=[dvec; col(d)]; %#ok<AGROW>
    corrvec=[corrvec; col(Corr)]; %#ok<AGROW>
    
    if length(w)>1 || ~isempty(problem)
        H2vec(i)=Inf;
        disconnectedcount=disconnectedcount+numsigs_permove;
    else
        
        sigmay=lyap(Ltilde,Q*noise*transpose(Q));
        H=sqrt(trace(sigmay));
        H2vec(i)=H;
    end
    
    receivers=randsample(N,numsigs_permove,'true');
    for j=1:numsigs_permove        

        receiver=receivers(j);
        allreceivers=d(receiver,:)<=radius;
        beta=zeros(N,1);
        beta(allreceivers)=b;
        B=diag(beta);
        Lf=L-B;
        Pf=P-B;
        [~,vals]=eig(Pf);
        problem2=find(sigfig(diag(vals),13)==0,1);
        Lftilde=Qfull*Lf*transpose(Qfull);
        Bbar=Qfull*B;
        [~,vals]=eig(Lf);
        problem3=find(sigfig(diag(vals),13)==0,1);

        if ~isempty(problem3) || ~isempty(problem2)
            nooneislistening=nooneislistening+1;
            H2vec_forced(i,j)=Inf;
        else
            sigmaz_forced=lyap(Lftilde,Qfull*noise*transpose(Qfull));
            sigmay_forced=sigmaz_forced(1:end-1,1:end-1);
            H_forced=sqrt(trace(sigmay_forced));
%             if imag(H_forced)==0
            H2vec_forced(i,j)=H_forced;
%             else H2vec_forced(i,j)=Inf;
%             end
        end
        
            toinvert=-Pf-1/(b*sum(allreceivers))*B*ones(N,1)*ones(1,N)*B;
            [W,Lambda]=eig(toinvert);
            w3=find(sigfig(diag(Lambda),13)==0);
            Cov_forced=W(:,setdiff(1:N,w3))*inv(Lambda(setdiff(1:N,w3),setdiff(1:N,w3)))*transpose(W(:,setdiff(1:N,w3)));
            todivide=repmat(diag(Cov_forced),1,N).*repmat(reshape(diag(Cov_forced),1,[]),N,1);
            todivide(sigfig(Cov_forced,13)==0)=1;
            todivide=power(todivide,.5);
            Corr_forced=Cov_forced./todivide;
            dvec_forced=[dvec_forced; col(d)]; %#ok<AGROW>
            corrvec_forced=[corrvec_forced; col(Corr_forced)]; %#ok<AGROW>
    end
end

if disconnectedcount<nummoves
    H2vec=col(H2vec);
    meanH2=mean(H2vec(~isinf(H2vec)));
    H2vec_forced=col(H2vec_forced);
    meanH2_forced=mean(H2vec_forced(~isinf(H2vec_forced)));
else
    meanH2=Inf;
    meanH2_forced=Inf;
end
meanrho=mean(1./H2vec);
meanrho_forced=mean(1./H2vec_forced);

deltad=.05;
dbins=(0:deltad:1.4)+deltad/2;
Nbins=length(dbins);
avgcorrs=zeros(1,Nbins);
avgcorrs_forced=zeros(1,Nbins);

for i=1:Nbins
    w1=find(dvec>=dbins(i)-deltad/2);
    w2=find(dvec<dbins(i)+deltad/2);
    w=intersect(w1,w2);
    avgcorrs(i)=mean(corrvec(w));
end

for i=1:Nbins
    w1=find(dvec_forced>=dbins(i)-deltad/2);
    w2=find(dvec_forced<dbins(i)+deltad/2);
    w=intersect(w1,w2);
    avgcorrs_forced(i)=mean(corrvec_forced(w));
end

f=find(avgcorrs<0,1,'first');
if avgcorrs(1)==0
    corrlength=0;
elseif isempty(f)
    corrlength=max(dvec);
else
    corrlength=interp1(avgcorrs([f-1 f]),dbins([f-1 f]),0);
end

f=find(avgcorrs_forced<0,1,'first');
if avgcorrs_forced(1)==0
    corrlength=0;
elseif isempty(f)
    corrlength_forced=max(dvec);
else
    corrlength_forced=interp1(avgcorrs_forced([f-1 f]),dbins([f-1 f]),0);
end



