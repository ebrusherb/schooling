
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
%%

% function [meanH2, meanH2_forced, corrlength, corrlength_forced, disconnectedcount]=groupprops(strategy,numsigs_permove,nummoves,radius,b,T)

% N=max(size(strategy));
% close all
dvec=[];
corrvec=[];
dvec_forced=[];
corrvec_forced=[];
disconnectedcount=0;
H2vec=zeros(nummoves,1);
H2vec_forced=zeros(nummoves,numsigs_permove);

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
    Lbar=Qfull*L*transpose(Qfull);
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
    else
        
        sigma=lyap(Ltilde,Q*noise*transpose(Q));
%         sigmaz=lyap(Lbar,Qfull*noise*transpose(Qfull));
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
    
    
%         receivers=randsample(N,numsigs_permove,'true');
        for j=1:numsigs_permove        

%             receiver=receivers(j);
            receiver=1;
            allreceivers=d(receiver,:)<=radius;
            beta=zeros(N,1);
            beta(allreceivers)=b;
            B=diag(beta);
%             Lf=L-B;
            Lf=P-B;
            Lfbar=Q*Lf*transpose(Q);
            Pf=P-B; 
            Pftilde=-Qfull*Pf*transpose(Qfull);
            Lftilde=Qfull*Lf*transpose(Qfull);
            Bti=Qfull*B;
            [~,vals]=eig(Lf);
            w=find(sigfig(diag(vals),13)==0,1);
            if ~isempty(w)
                disconnectedcount=disconnectedcount+1;
                H2vec_forced(i,j)=H;
                dvec_forced=[dvec_forced; col(d)]; %#ok<AGROW>
                corrvec_forced=[corrvec_forced; col(C)];%#ok<AGROW>
            else
                sigma_forced=lyap(Lf,-inv(Lf)*B*ones(N,1)*ones(1,N)*B-B*ones(N,1)*ones(1,N)*B*transpose(inv(Lf))+noise);
                sigmaz=lyap(Lftilde,Qfull*transpose(Qfull));
                sigmaz2=Qfull*sigma_forced*transpose(Qfull);
                sigmaz3=lyap(Lftilde,-inv(Lftilde)*Bbar*ones(N,1)*ones(1,N)*transpose(Bbar)-Bbar*ones(N,1)*ones(1,N)*transpose(Bbar)*transpose(inv(Lftilde))+Qfull*transpose(Qfull));
                sigmay=sigmaz2(1:end-1,1:end-1);
                H_forced=sqrt(trace(sigmay));
%                 if trace(sigmay)>1e3
%                     'big H'
%                     figure
%                     subplot(2,2,1);
%                     imagesc(sigma);
% %                     imagesc(sigmaz2);
%                     colorbar;
%                     subplot(2,2,2);
%                     imagesc(sigmay);
% %                     imagesc(sigmaz3);
%                     colorbar;
%                     subplot(2,2,3);
%                     imagesc(sigma_forced);
%                     colorbar;
%                     subplot(2,2,4);
%                     imagesc(transpose(Qfull)*sigmaz3*Qfull);
%                     colorbar
%                     i=nummoves+1;
%                 end
                tocheck=Lftilde*sigmaz3+sigmaz3*transpose(Lftilde)-inv(Lftilde)*Bbar*ones(N,1)*ones(1,N)*transpose(Bbar)-Bbar*ones(N,1)*ones(1,N)*transpose(Bbar)*transpose(inv(Lftilde))+Qfull*transpose(Qfull);
                m=max(max(abs(tocheck(1:end-1,1:end-1))));
                if m<1e-12
                    H2vec_forced(i,j)=H_forced;

    %                 R=-Qfull*Pf*transpose(Qfull);
    %                 toinvert=-Pf-1/R(1,1)*Pf*1/N*ones(N,N)*Pf;
    %                 R11=b/N*sum(allreceivers);
    %                 toinvert=-Pf-1/R11*Pf*1/N*ones(N,N)*Pf;
                    toinvertz=-Pf-1/(b*sum(allreceivers))*B*ones(N,1)*ones(1,N)*B;
                    toinvert_y=Pftilde-N/(b*sum(allreceivers))*Pftilde(:,end)*Pftilde(end,:);
                    toinvert_y=toinvert_y(1:end-1,1:end-1);
                    [vecs,vals]=eig(toinvertz);
                    w=find(sigfig(diag(vals),13)==0);
                    inverted=vecs(:,[1:(w-1) (w+1):end])*inv(vals([1:(w-1) (w+1):end],[1:(w-1) (w+1):end]))*transpose(vecs(:,[1:(w-1) (w+1):end]));
                    todivide=repmat(diag(inverted),1,N).*repmat(reshape(diag(inverted),1,[]),N,1);
                    todivide=power(todivide,.5);
                    Cforced=inverted./todivide;
                    dvec_forced=[dvec_forced; col(d)]; %#ok<AGROW>
                    corrvec_forced=[corrvec_forced; col(Cforced)]; %#ok<AGROW>
                else
                    disconnectedcount=disconnectedcount+1;
                    H2vec_forced(i,j)=H;
                    dvec_forced=[dvec_forced; col(d)]; %#ok<AGROW>
                    corrvec_forced=[corrvec_forced; col(C)];%#ok<AGROW>
                end
            end
        end
    end
    i=i+1;
end

meanH2=mean(H2vec(~isnan(H2vec)));
H2vec_forced=col(H2vec_forced);
meanH2_forced=mean(H2vec_forced(~isnan(H2vec_forced)));
% meanH2_forced=meanH2;

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

% [~,vals]=eig(Lftilde);
% % corrlength_forced=corrlength;
% [S,~]=graphconncomp(sparse(A),'Directed','true');
% [S min(abs(diag(vals)))]
% % 
% % [S,~]=graphconncomp(sparse(A),'Directed','true');
% % figure
% % subplot(1,2,1)
% % imagesc(Ltilde*sigma+sigma*transpose(Ltilde)+Q*transpose(Q))
% % colorbar
% % subplot(1,2,2)
% % toplot=Lftilde*sigmaz3+sigmaz3*transpose(Lftilde)-inv(Lftilde)*Bbar*ones(N,1)*ones(1,N)*transpose(Bbar)-Bbar*ones(N,1)*ones(1,N)*transpose(Bbar)*transpose(inv(Lftilde))+Qfull*transpose(Qfull);
% % imagesc(toplot(1:end-1,1:end-1));
% % % imagesc(Lftilde*sigmaz3+sigmaz3*transpose(Lftilde)-inv(Lftilde)*Bbar*ones(N,1)*ones(1,N)*transpose(Bbar)-Bbar*ones(N,1)*ones(1,N)*transpose(Bbar)*transpose(inv(Lftilde))+Qfull*transpose(Qfull));
% % colorbar

subplot(1,3,1);
% imagesc(sigmay);colorbar
imagesc(transpose(Q)*sigmay*Q);colorbar;
% imagesc(Lfbar*sigmay+sigmay*transpose(Lfbar));colorbar
% imagesc(sigmaz2);colorbar
subplot(1,3,2);
% imagesc(inv(toinvert_y));colorbar
% imagesc(transpose(Q)*inv(toinvert_y)*Q);colorbar
% imagesc(Lfbar*inv(toinvert_y)+inv(toinvert_y)*Lfbar);colorbar
imagesc(inverted);colorbar
subplot(1,3,3);
% plot(sigmay,inv(toinvert_y),'o')
plot(transpose(Q)*sigmay*Q,inverted,'o')
%%
Pftilde=-Qfull*Pf*transpose(Qfull);
inv_zcov=-Pf-1/(b*sum(allreceivers))*B*ones(N,1)*ones(1,N)*B;
inv_ycov=Pftilde-1/Pftilde(N,N)*Pftilde(:,end)*Pftilde(end,:);
[vecs,vals]=eig(inv_zcov);
w=find(sigfig(diag(vals),13)==0);
vals=diag(vals);
vecs=vecs(:,[1:(w-1) (w+1):end]);
vals=vals([1:(w-1) (w+1):end]);

zcov=vecs*inv(diag(vals))*transpose(vecs);

[det(-inv(Pf))*Pftilde(N,N) det(inv(inv_ycov(1:end-1,1:end-1))) 1/prod(vals)]


%%
eN=[zeros(19,1);1];
v=rand(N,1);
V=mean(v);
z=v-V;
y=Qfull*z;
yprime=y+V*sqrt(N)*eN;
[v transpose(Qfull)*yprime]

%%
Pftilde=-Qfull*Pf*transpose(Qfull);
v=rand(N,1);
V=mean(v);
z=v-V;
y=Qfull*z;
pv=probv(v,Pf)
pzV=probzV(z,V,Pf)
py=probyV(y,V,Pf,Pftilde)

xvals=-10:.1:10;
yvals=-10:.1:10;
L=length(xvals);
probmat_v=zeros(L,L);
probmat_z=zeros(L,L);
probmat_y=zeros(L,L);
probvec_y=zeros(1,L);
for i=1:L
    Y=xvals(i);
    probvec_y(i)=proby([Y;0],Pf,Pftilde);
    for j=1:L
        V2=yvals(j);
        y=[Y;0];
        v=[yvals(i);yvals(j)];
        V=mean(v);
        z=v-V;
        probmat_v(i,j)=probv(v,Pf);
        probmat_z(i,j)=probzV(z,V,Pf);
        probmat_y(i,j)=probyV(y,V2,Pf,Pftilde);
    end
end


subplot(1,3,1)
imagesc(xvals,yvals,probmat_v);set(gca,'ydir','normal');colorbar
subplot(1,3,2)
imagesc(xvals,yvals,probmat_z);set(gca,'ydir','normal');colorbar
subplot(1,3,3)
imagesc(xvals,yvals,probmat_y);set(gca,'ydir','normal');colorbar