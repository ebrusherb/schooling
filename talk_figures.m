figure
set(gcf,'Color','w')
w=6.83;
h=.5*w;
set(gcf,'Units','inches');
set(gcf,'Position',[11.5 3 w h]);

N=20;
strategy=3*ones(N,1);
positions=unifrnd(0,1,N,2);
d=squareform(pdist(positions));

M=zeros(N);
for ind=1:N
    [~, order]=sort(d(ind,:));
    neighbors=order(2:strategy(ind)+1);
    M(ind,neighbors)=1/strategy(ind);
    M(ind,ind)=-sum(M(ind,neighbors));
end
tograph=sparse(M);
gplot(tograph,positions,'-o')
f=findobj(gcf,'Marker','o');
set(f,'MarkerSize',8,'Color','w','MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',3)

hold on
delta=0.8;
delta2=0.2;
arrowlength=10;arrow_intangle=15;arrow_angle=30;lw=2;
for i=1:N
    for j=setdiff(1:N,i)
        if M(i,j)~=0
            v=positions(j,:)-positions(i,:);
            start=positions(j,:)-delta*v;
            arrowend=positions(j,:)-delta2*v;
            arrow(start,arrowend,5,15,30,'EdgeColor','w','MarkerEdgeColor','w','LineWidth',lw)
        end
    end
end
box off
gplot(tograph,positions,'-o')
f=findobj(gcf,'Marker','o');
set(f,'MarkerSize',8,'Color','w','MarkerFaceColor','r','MarkerEdgeColor','r','LineWidth',3)

axis off

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename='/Users/eleanorbrush/Desktop/network6.pdf';
print(filename,'-dpdf','-r300')

%%
N=20;radius=.1;b=1;T=1;nummoves=1000;numsigs_permove=1;
strategy=6*ones(N,1);


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

positions=unifrnd(0,1,N,2);
    d=squareform(pdist(positions));
    [~,o]=sort(d(1,:));
    d=d(o,o);

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
    noise=diag(strategy);
    
    S=.5*(Abar+transpose(Abar));
    P=S;
    P(1:N+1:end)=-sum(S,2);
    
    receiver=1;
    allreceivers=d(receiver,:)<=radius;
    beta=zeros(N,1);
    beta(allreceivers)=b;
    B=diag(beta);
    Pf=P-B; 

    toinvert=-Pf-1/(b*sum(allreceivers))*B*ones(N,1)*ones(1,N)*B;
    [W,Lambda]=eig(toinvert);
    w3=find(sigfig(diag(Lambda),13)==0);
    Cov_forced=W(:,setdiff(1:N,w3))*inv(Lambda(setdiff(1:N,w3),setdiff(1:N,w3)))*transpose(W(:,setdiff(1:N,w3)));
    todivide=repmat(diag(Cov_forced),1,N).*repmat(reshape(diag(Cov_forced),1,[]),N,1);
    todivide=power(todivide,.5);
    Corr_forced=Cov_forced./todivide;
    
figure
set(gcf,'Color','w')
w=4;
h=w;
set(gcf,'Units','inches');
set(gcf,'Position',[11.5 3 w h]);

divcols=cbrewer('div','RdBu',20);
imagesc(-Corr_forced);
colormap(divcols);
caxis manual 
caxis([-1 1]);


filename='/Users/eleanorbrush/Desktop/corr_heatmap.pdf';
print(filename,'-dpdf','-r300')

    
