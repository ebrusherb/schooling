
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
close all
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
    S=.5*(Abar+transpose(Abar));
    P=S;
    P(1:N+1:end)=-sum(S,2);
    L=P;
    
    Lbar=Q*L*transpose(Q);
    Ltilde=Qfull*L*transpose(Qfull);
    Pbar=Q*P*transpose(Q);
    
    [~,vals]=eig(Lbar);
    problem=find(sigfig(diag(vals),13)==0,1);
    
    noise=eye(N);
%     noise=diag(strategy);
    
    [W,Lambda]=eig(P);
    w=find(sigfig(diag(Lambda),13)==0);
    
    if length(w)>1 || ~isempty(problem)
        disconnectedcount=disconnectedcount+numsigs_permove;
    else
        
        sigmay=lyap(Lbar,Q*noise*transpose(Q));
%         sigmaz=lyap(Lbar,Qfull*noise*transpose(Qfull));
        H=sqrt(trace(sigmay));
        H2vec(i)=H;
      
        Wtilde=W(:,[1:(w-1) (w+1):end]);
        Lambdatilde=Lambda([1:(w-1) (w+1):end],[1:(w-1) (w+1):end]);
        Cov=Wtilde*inv(Lambdatilde)*transpose(Wtilde); %#ok<*MINV>
        Cy=-Q*Cov*transpose(Q);
        
        todivide=repmat(diag(-Cov),1,N).*repmat(reshape(diag(-Cov),1,[]),N,1);
        todivide=power(todivide,.5);
        Corr=-Cov./todivide;
        dvec=[dvec; col(d)]; %#ok<AGROW>
        corrvec=[corrvec; col(Corr)]; %#ok<AGROW>
    
    
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
            Lftilde=Qfull*Lf*transpose(Qfull);
            Pf=P-B; 
            Pfbar=Q*Pf*transpose(Q);
            Pftilde=-Qfull*Pf*transpose(Qfull);
            Bbar=Qfull*B;
            [~,vals]=eig(Lf);
            w=find(sigfig(diag(vals),13)==0,1);
            if ~isempty(w)
                disconnectedcount=disconnectedcount+1;
                H2vec_forced(i,j)=H;
                dvec_forced=[dvec_forced; col(d)]; %#ok<AGROW>
                corrvec_forced=[corrvec_forced; col(Corr)];%#ok<AGROW>
            else
                sigma_ybar3=lyap(Lf,-inv(Lf)*B*ones(N,1)*ones(1,N)*B-B*ones(N,1)*ones(1,N)*B*transpose(inv(Lf))+noise);
                sigma_ybar=lyap(Lftilde,Qfull*noise*transpose(Qfull));
%                 sigma_ybar2=Qfull*sigma_ybar3*transpose(Qfull);
%                 sigma_ybar4=lyap(Lftilde,-inv(Lftilde)*Bbar*ones(N,1)*ones(1,N)*transpose(Bbar)-Bbar*ones(N,1)*ones(1,N)*transpose(Bbar)*transpose(inv(Lftilde))+Qfull*transpose(Qfull));
                sigmay_forced=sigma_ybar(1:end-1,1:end-1);
                H_forced=sqrt(trace(sigmay_forced));

                H2vec_forced(i,j)=H_forced;

                toinvert_z=-Pf-1/(b*sum(allreceivers))*B*ones(N,1)*ones(1,N)*B;
                toinvert_y=Pftilde-N/(b*sum(allreceivers))*Pftilde(:,end)*Pftilde(end,:);
                toinvert_y=toinvert_y(1:end-1,1:end-1);
                Cy_forced=inv(toinvert_y);
%                 Cy_forced=Q*inv(-Pf-1/(b*sum(allreceivers))*B*ones(N,1)*ones(1,N)*B)*transpose(Q);
                [W,Lambda]=eig(toinvert_z);
                w=find(sigfig(diag(Lambda),13)==0);
                inverted_z=W(:,[1:(w-1) (w+1):end])*inv(Lambda([1:(w-1) (w+1):end],[1:(w-1) (w+1):end]))*transpose(W(:,[1:(w-1) (w+1):end]));
                Cov_forced=inverted_z;
                todivide=repmat(diag(Cov_forced),1,N).*repmat(reshape(diag(Cov_forced),1,[]),N,1);
                todivide=power(todivide,.5);
                Corr_forced=Cov_forced./todivide;
                dvec_forced=[dvec_forced; col(d)]; %#ok<AGROW>
                corrvec_forced=[corrvec_forced; col(Corr_forced)]; %#ok<AGROW>

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

figure
subplot(2,2,1);
imagesc(sigmay);
% imagesc(Lbar*sigmay+sigmay*transpose(Lbar));
colorbar
subplot(2,2,2);
imagesc(Cy);
% imagesc(Pbar*Cy+Cy*transpose(Pbar));
colorbar
subplot(2,2,3);
imagesc(sigmay_forced);
% imagesc(Lfbar*sigmay_forced+sigmay_forced*transpose(Lfbar));
colorbar
subplot(2,2,4);
imagesc(Cy_forced);
% imagesc(Pfbar*Cy_forced+Cy_forced*transpose(Pfbar));
colorbar

figure
subplot(1,2,1)
plot(sigmay,Cy,'o')
subplot(1,2,2)
plot(sigmay_forced,Cy_forced,'o')
%%
subplot(1,2,1)
imagesc(sigma_ybar);
colorbar
subplot(1,2,2)
imagesc(Cy_forced)
colorbar
%%

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

%%
m1=transpose(Qfull)*(Pftilde-1/Pftilde(N,N)*Pftilde(:,end)*Pftilde(end,:))*Qfull;
m2=-Pf-1/(b*sum(allreceivers))*B*ones(N,1)*ones(1,N)*B;
% m1=transpose(Qfull)*Pftilde*Qfull;
% m2=-Pf;
