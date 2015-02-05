% function [meanH2, meanH2_forced, corrlength, corrlength_forced, disconnectedcount]=groupprops(strategy,numsigs_permove,nummoves,radius,b,T)

N=max(size(strategy));

dvec=[];
corrvec=[];
dvec_forced=[];
corrvec_forced=[];
disconnectedcount=0;
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
        H=trace(sigma);
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
        for j=1:numsigs_permove        

            receiver=receivers(j);
            allreceivers=d(receiver,:)<=radius;
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
            else
                sigma_forced=lyap(Lf,eye(N));
                H_forced=trace(sigma_forced);
                H2vec_forced(i,j)=H_forced;

                R=-Qfull*Pf*transpose(Qfull);
                toinvert=-Pf-1/R(1,1)*Pf*1/N*ones(N,N)*Pf;
                [vecs,vals]=eig(toinvert);
                w=find(sigfig(diag(vals),13)==0);
                inverted=vecs(:,[1:(w-1) (w+1):end])*inv(vals([1:(w-1) (w+1):end],[1:(w-1) (w+1):end]))*transpose(vecs(:,[1:(w-1) (w+1):end]));
                todivide=repmat(diag(inverted),1,N).*repmat(reshape(diag(inverted),1,[]),N,1);
                todivide=power(todivide,.5);
                C=inverted./todivide;
                dvec_forced=[dvec_forced; col(d)]; %#ok<AGROW>
                corrvec_forced=[corrvec_forced; col(C)]; %#ok<AGROW>
            end
        end
    end
end

meanH2=mean(H2vec(~isnan(H2vec)));
H2vec_forced=col(H2vec_forced);
meanH2_forced=mean(H2vec_forced(~isnan(H2vec_forced)));

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


