N=20;
b=1;
T=1;
nummoves=5000;
numsigs_permove=1;

strategy=6*ones(N,1);
radius=1.1;

storeP=zeros(nummoves,N,N);
storereceivers=zeros(nummoves,N);
dvec=[];
corrvec=[];
dvec_forced=[];
corrvec_forced=[];
disconnectedcount=0;
H2vec=zeros(nummoves,1);
H2vec_forced=zeros(nummoves);

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
Qfull=[1/sqrt(N)*ones(1,N); Q];

i=1;
while i<=nummoves   
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
    
    noise=eye(N);
%     noise=diag(strategy);
    
    S=.5*(Abar+transpose(Abar));
    P=S;
    P(1:N+1:end)=-sum(S,2);
    [W,Lambda]=eig(P);
    w=find(sigfig(diag(Lambda),13)==0);
    
    
    
    if length(w)>1 || ~isempty(problem)
        disconnectedcount=disconnectedcount+numsigs_permove;
        dvec=[dvec; col(d)]; %#ok<AGROW>
        corrvec=[corrvec; zeros(N^2,1)]; %#ok<AGROW>
        dvec_forced=[dvec_forced; col(d)]; %#ok<AGROW>
        corrvec_forced=[corrvec_forced; zeros(N^2,1)]; %#ok<AGROW>
        i=i+1;
    else
        storeP(i,:,:)=P;
        sigma=lyap(Ltilde,Q*noise*transpose(Q));
        H=sqrt(trace(sigma));
        H2vec(i)=H;
        
%         sigma_forced=lyap(Lf,eye(N));
%         H_forced=trace(sigma_forced);
%         H2vec_forced(i)=H_forced;
%         
        Wtilde=W(:,[1:(w-1) (w+1):end]);
        Lambdatilde=Lambda([1:(w-1) (w+1):end],[1:(w-1) (w+1):end]);
        Ptildeinv=Wtilde*inv(Lambdatilde)*transpose(Wtilde); %#ok<*MINV>
        
        todivide=repmat(diag(-Ptildeinv),1,N).*repmat(reshape(diag(-Ptildeinv),1,[]),N,1);
        todivide=power(todivide,.5);
        C=-Ptildeinv./todivide;
        dvec=[dvec; col(d)]; %#ok<AGROW>
        corrvec=[corrvec; col(C)]; %#ok<AGROW>
    
    
        receivers=randsample(N,numsigs_permove,'true');
        

            receiver=receivers(1);
            allreceivers=d(receiver,:)<=radius;
            storereceivers(i,:)=allreceivers;
            
            beta=zeros(N,1);
            beta(allreceivers)=b;
            B=diag(beta);
            Lf=L-B;
            Pf=P-B; 
            [~,vals]=eig(Lf);
            w=find(sigfig(diag(vals),13)==0,1);
            if ~isempty(w)
                disconnectedcount=disconnectedcount+1;
                H2vec_forced(i,j)=H;
                dvec_forced=[dvec_forced; col(d)]; %#ok<AGROW>
                corrvec_forced=[corrvec_forced; corrvec];%#ok<AGROW>
                i=i+1;
            else
                sigma_forced=lyap(Lf,noise);
                H_forced=trace(sigma_forced);
                H2vec_forced(i)=H_forced;

%                 R=-Qfull*Pf*transpose(Qfull);
%                 toinvert=-Pf-1/R(1,1)*Pf*1/N*ones(N,N)*Pf;
%                 R11=b/N*sum(allreceivers);
%                 toinvert=-Pf-1/R11*Pf*1/N*ones(N,N)*Pf;
                toinvert=-Pf-1/(b*sum(allreceivers))*B*ones(N,1)*ones(1,N)*B;
                [vecs,vals]=eig(toinvert);
                w=find(sigfig(diag(vals),13)==0);
                inverted=vecs(:,[1:(w-1) (w+1):end])*inv(vals([1:(w-1) (w+1):end],[1:(w-1) (w+1):end]))*transpose(vecs(:,[1:(w-1) (w+1):end]));
                todivide=repmat(diag(inverted),1,N).*repmat(reshape(diag(inverted),1,[]),N,1);
                todivide=power(todivide,.5);
                C=inverted./todivide;
                dvec_forced=[dvec_forced; col(d)]; %#ok<AGROW>
                corrvec_forced=[corrvec_forced; col(C)]; %#ok<AGROW>
                if max(max(abs(imag(C))))~=0
                    i %#ok<NOPTS>
                    'theres a problem' %#ok<NOPTS>
                    i=nummoves+1;
                else i=i+1;
                end
            end
        
    end
end

meanH2=mean(H2vec(~isnan(H2vec)));
H2vec_forced=col(H2vec_forced);
meanH2_forced=mean(H2vec_forced(~isnan(H2vec_forced)));
% meanH2_forced=meanH2;
%%
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
if isempty(f)
    corrlength=max(dvec);
else
    corrlength=interp1(avgcorrs([f-1 f]),dbins([f-1 f]),0);
end

f=find(avgcorrs_forced<0,1,'first');
if isempty(f)
    corrlength_forced=max(dvec);
else
    corrlength_forced=interp1(avgcorrs_forced([f-1 f]),dbins([f-1 f]),0);
end

corrlength_forced=corrlength;
%%
load('/Users/eleanorbrush/Desktop/problematic.mat')
R=-Qfull*Pf*transpose(Qfull);
toinvert=-Pf-1/R(1,1)*Pf*1/N*ones(N,N)*Pf;
[vecs,vals]=eig(toinvert);
R11=b/N*sum(allreceivers);
toinvert2=-Pf-1/R11*Pf*1/N*ones(N,N)*Pf;
[vecs2,vals2]=eig(toinvert2);
% vals=sigfig(vals,14);
w=find(sigfig(diag(vals),13)==0);
inverted=vecs(:,[1:(w-1) (w+1):end])*inv(vals([1:(w-1) (w+1):end],[1:(w-1) (w+1):end]))*transpose(vecs(:,[1:(w-1) (w+1):end]));
todivide=repmat(diag(inverted),1,N).*repmat(reshape(diag(inverted),1,[]),N,1);
todivide=power(todivide,.5);
C=inverted./todivide;

%%
f=find(imag(corrvec_forced)~=0);
[j,k,i]=ind2sub([N,N,nummoves],f);
ivals=unique(i);
%%
i=ivals(2);
P=reshape(storeP(i,:,:),N,N);
allreceivers=storereceivers(i,1,:);
beta=zeros(N,1);
beta(allreceivers)=b;
B=diag(beta);
Pf=P-B;
%%
R=-Qfull*Pf*transpose(Qfull);
toinvert=-Pf-1/R(1,1)*Pf*1/N*ones(N,N)*Pf;
[vecs,vals]=eig(toinvert);
w=find(sigfig(diag(vals),13)==0);
inverted=vecs(:,[1:(w-1) (w+1):end])*inv(vals([1:(w-1) (w+1):end],[1:(w-1) (w+1):end]))*transpose(vecs(:,[1:(w-1) (w+1):end]));
todivide=repmat(diag(inverted),1,N).*repmat(reshape(diag(inverted),1,[]),N,1);
todivide=power(todivide,.5);
C=inverted./todivide;
