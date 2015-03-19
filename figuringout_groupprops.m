strategy=3*ones(N,1);

N=max(size(strategy));

dvec=[];
corrvec=[];
corrvec_old=[];
% dvec_forced=[];
% corrvec_forced=[];
disconnectedcount=0;
H2vec=zeros(nummoves,1);
H2vec_old=zeros(nummoves,1);
% H2vec_forced=zeros(nummoves,numsigs_permove);

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
    
    S=.5*(Abar+transpose(Abar));
    P=S;
    P(1:N+1:end)=-sum(S,2);
    [W,Lambda]=eig(P);
    w=find(sigfig(diag(Lambda),13)==0);
    
    if length(w)>1 || ~isempty(problem)
        disconnectedcount=disconnectedcount+numsigs_permove;
    else
        
        sigma=lyap(Ltilde,eye(N-1));
        H=sqrt(trace(sigma));
        H2vec(i)=H;
        
        h=H2norm(L,'non');
        H2vec_old(i)=h; 
        
        Wtilde=W(:,[1:(w-1) (w+1):end]);
        Lambdatilde=Lambda([1:(w-1) (w+1):end],[1:(w-1) (w+1):end]);
        Ptildeinv=Wtilde*inv(Lambdatilde)*transpose(Wtilde); %#ok<*MINV>
        
        todivide=repmat(diag(-Ptildeinv),1,N).*repmat(reshape(diag(-Ptildeinv),1,[]),N,1);
        todivide=power(todivide,.5);
%         C=-Ptildeinv./todivide;
        C=-Ptildeinv;
        dvec=[dvec; col(d)]; %#ok<AGROW>
        corrvec=[corrvec; col(C)]; %#ok<AGROW>
        
        corrs=paircorrelations(L,zeros(N,1));
        corrvec_old=[corrvec_old; col(corrs)]; %#ok<AGROW>
    
    end
end

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