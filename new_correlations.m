N=20;
b=2;
informed=15;
radius=.5;
strategy=[10,10*ones(1,N-1)];

dvec=[];
corrvec=[];
corrvec_forced=[];
disconnectedcount=0;
nummoves=100;
H2vec1=zeros(1,nummoves);
H2vec_forced1=zeros(1,nummoves);

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
    prob=find(sigfig(diag(vals),13)==0);
    
    S=.5*(Abar+transpose(Abar));
    P=S;
    P(1:N+1:end)=-sum(S,2);
    [W,Lambda]=eig(P);
    w=find(sigfig(diag(Lambda),13)==0);
    
    if length(w)>1 || ~isempty(prob)
        disconnectedcount=disconnectedcount+1;
    else
        
        sigma=lyap(Ltilde,eye(N-1));
        H=trace(sigma);
        H2vec1(i)=H;
         
        Wtilde=W(:,[1:(w-1) (w+1):end]);
        Lambdatilde=Lambda([1:(w-1) (w+1):end],[1:(w-1) (w+1):end]);
        Ptildeinv=Wtilde*inv(Lambdatilde)*transpose(Wtilde);
        
        todivide=repmat(diag(-Ptildeinv),1,N).*repmat(reshape(diag(-Ptildeinv),1,[]),N,1);
        todivide=power(todivide,.5);
        C=-Ptildeinv./todivide;
        dvec=[dvec; col(d)];
        corrvec=[corrvec; col(C)];
        
    end
   
    
%     receiver=randi(20,1);
    receiver=randsample(20,1);
    allreceivers=d(receiver,:)<=radius;
    beta=zeros(N,1);
    beta(allreceivers)=b;
    B=diag(beta);
    
    Lf=L-B;
    Pf=P-B;

    sigma_forced=lyap(Lf,eye(N));
    H_forced=trace(sigma_forced);
    H2vec_forced1(i)=H_forced;

    R=-Qfull*Pf*transpose(Qfull);
    toinvert=-Pf-1/R(1,1)*Pf*1/N*ones(N,N)*Pf;
    [vecs,vals]=eig(toinvert);
    w=find(sigfig(diag(vals),13)==0);
    inverted=vecs(:,[1:(w-1) (w+1):end])*inv(vals([1:(w-1) (w+1):end],[1:(w-1) (w+1):end]))*transpose(vecs(:,[1:(w-1) (w+1):end]));
    todivide=repmat(diag(inverted),1,N).*repmat(reshape(diag(inverted),1,[]),N,1);
    todivide=power(todivide,.5);
    C=inverted./todivide;
    corrvec_forced=[corrvec_forced; col(C)];
end

meanH2=mean(H2vec1(~isnan(H2vec1)));
meanH2_forced=mean(H2vec_forced1(~isnan(H2vec_forced1)));


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



%%

R=-Qfull*Pf*transpose(Qfull);
check1=inv(Qfull)*(R-1/R(1,1)*R(:,1)*R(1,:))*transpose(inv(Qfull));
hold1=Pf*transpose(Qfull);
hold2=Qfull*Pf;
% check2=-Pf-1/R(1,1)*hold1(:,1)*hold2(1,:);
check2=-Pf-1/R(1,1)*Pf*transpose(Qfull(1,:))*Qfull(1,:)*Pf;
subplot(1,3,1)
imagesc(check1)
colorbar
subplot(1,3,2)
imagesc(check2)
colorbar
toinvert=-Pf-1/R(1,1)*Pf*1/N*ones(N,N)*Pf;
subplot(1,3,3)
[vecs,vals]=eig(toinvert);
w=find(sigfig(diag(vals),13)==0);
if length(w)==1
    inverted=vecs(:,[1:(w-1) (w+1):end])*inv(vals([1:(w-1) (w+1):end],[1:(w-1) (w+1):end]))*transpose(vecs(:,[1:(w-1) (w+1):end]));
end
%%
sqrt(2*pi)*sqrt(R(1,1)*N)/power(sqrt(2*pi),N/2)/sqrt(det(-inv(Pf)))*power(2*pi,N/2)%*det(inverted)
%%
Qfull=[1/sqrt(N)*ones(1,N); Q];
R=-Qfull*Pf*transpose(Qfull);
v=rand(N,1);
V=mean(v);
z=v-V;
y=transpose(inv(Qfull))*z;
e1=[1;zeros(N-1,1)];
exp(-1/2*transpose(v-ones(N,1))*-Pf*(v-ones(N,1)))
exp(-1/2*transpose(z-(1-V)*ones(N,1))*-Pf*(z-(1-V)*ones(N,1)))
exp(-1/2*transpose(y-(1-V)*sqrt(N)*e1)*R*(y-(1-V)*sqrt(N)*e1))
exp(-1/2*(transpose(y)*R*y-power(R(1,:)*y/sqrt(R(1,1)),2)))*exp(-1/2*(R(1,1)*N*(V-1)^2+2*(V-1)*sqrt(N)*R(1,:)*y+power(R(1,:)*y/sqrt(R(1,1)),2)))
exp(-1/2*(transpose(y)*R*y-power(R(1,:)*y/sqrt(R(1,1)),2)))*exp(-1/2*(sqrt(R(1,1)*N)*(V-1)+R(1,:)*y/sqrt(R(1,1)))^2)
exp(-1/2*(transpose(y)*(R-R(:,1)*R(1,:)/((R(1,1))))*y))*exp(-1/2*(sqrt(R(1,1)*N)*(V-1)+R(1,:)*y/sqrt(R(1,1)))^2)

%%
toinvert=-Pf-1/R(1,1)*Pf*1/N*ones(N,N)*Pf;
