
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
%     L=P;
    
    Lbar=Q*L*transpose(Q);
    Ltilde=Qfull*L*transpose(Qfull);
    Pbar=Q*P*transpose(Q);
    Ptilde=Qfull*P*transpose(Q);
    
    [~,vals]=eig(Lbar);
    problem=find(sigfig(diag(vals),13)==0,1);
    
    noise=eye(N);
%     noise=diag(strategy);
    
    [~,Lambda]=eig(P);
    w=find(sigfig(diag(Lambda),13)==0);
    
    if length(w)>1 || ~isempty(problem)
        disconnectedcount=disconnectedcount+numsigs_permove;
    else
        
        sigma_y=lyap(Lbar,Q*noise*transpose(Q));
        H=sqrt(trace(sigma_y));
        H2vec(i)=H;
      
        Pdagger=dagger(P);
        Cov=Pdagger; %#ok<*MINV>
        C_y=-Q*Cov*transpose(Q);
        
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
            Lf=L-B;
%             Lf=P-B;
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
                sigma_ybar=lyap(Lftilde,Qfull*noise*transpose(Qfull));
                sigma_y_forced=sigma_ybar(1:end-1,1:end-1);
                H_forced=sqrt(trace(sigma_y_forced));

                H2vec_forced(i,j)=H_forced;

                toinvert_z=-Pf-1/(b*sum(allreceivers))*B*ones(N,1)*ones(1,N)*B;
        
                Cov_forced=dagger(toinvert_z);
                
                C_y_forced=dagger(Pftilde-1/Pftilde(N,N)*Pftilde(:,end)*Pftilde(end,:));
                C_y_forced2=Q*Cov_forced*transpose(Q);
                
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
% imagesc(sigma_y);
imagesc(Lbar*sigma_y+sigma_y*transpose(Lbar));
colorbar
subplot(2,2,2);
% imagesc(C_y);
imagesc(Lbar*C_y+C_y*transpose(Lbar));
colorbar
subplot(2,2,3);
imagesc(sigma_y_forced);
% imagesc(Lfbar*sigma_y_forced+sigma_y_forced*transpose(Lfbar));
% imagesc(Lftilde*sigma_ybar+sigma_ybar*transpose(Lftilde));
colorbar
subplot(2,2,4);
imagesc(C_y_forced);
% imagesc(Pftilde*C_y_forced+C_y_forced(1:end-1,1:end-1)*transpose(Pfbar));
% imagesc(Lftilde*C_y_forced+C_y_forced*transpose(Lftilde));
colorbar

figure
subplot(1,2,1)
plot(sigma_y,C_y,'o')
hold on;plot(get(gca,'xlim'),2*get(gca,'xlim'))
subplot(1,2,2)
plot(sigma_y_forced,C_y_forced(1:end-1,1:end-1),'o')
% plot(sigma_ybar,C_y_forced,'o')
hold on;plot(get(gca,'xlim'),2*get(gca,'xlim'))
%%
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=.25*w;
set(gcf,'Units','inches');
set(gcf,'Position',[11.5 3 w h]);

[~,o]=sort(d(receiver,:));

xoffset=.08;
yoffset=.03;

subplot(1,3,1)
toplot=dagger(-Pf);
imagesc(toplot(o,o))
set(gca,'FontName',fontname,'FontSize',labfontsz,'tickdir','out')
axis square
cb=colorbar;
set(cb,'FontName',fontname,'FontSize',labfontsz)
apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String','(a)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

subplot(1,3,2)
imagesc(Cov_forced(o,o))
set(gca,'FontName',fontname,'FontSize',labfontsz,'tickdir','out')
axis square
cb=colorbar;
set(cb,'FontName',fontname,'FontSize',labfontsz)
apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String','(b)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

subplot(1,3,3)
imagesc(Corr_forced(o,o))
set(gca,'FontName',fontname,'FontSize',labfontsz,'tickdir','out')
axis square
cb=colorbar;
set(cb,'FontName',fontname,'FontSize',labfontsz)
apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String','(c)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);
filename=strcat('/Users/eleanorbrush/Desktop/','cov_matrices','.pdf');
print(filename,'-dpdf','-r300');
%%
%checking determinants make sense
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
%checking probabilities make sense
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

sum(sum(probmat_v))*.1*.1
sum(sum(probmat_z))*.1*.1
sum(sum(probmat_y))*.1*.1

