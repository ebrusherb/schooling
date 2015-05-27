N=20;
k2=17;
k1=6;
% radius=.1;
strategy=k2*ones(N-1,1);
strategy=[k1;strategy];
nummoves=1000;
m=2;
% p=(m)/N; 
freqs=zeros(N,nummoves);
ranks=zeros(N,nummoves);
for i=1:nummoves
%     draw=binornd(strategy,p);
    draw=hygernd(N,m,strategy);
    freq=draw./strategy;
    freqs(:,i)=freq;
    [~,~,rank]=unique(-freq);
    ranks(:,i)=rank;
end
eaten=zeros(N,nummoves);
gettoeat=zeros(N,nummoves);
for i=1:N
   eaten(i,:)=(ranks(i,:)==max(ranks));
   gettoeat(i,:)=(ranks(i,:)==min(ranks));
end
eaten=sum(eaten,2);
gettoeat=sum(gettoeat,2);
% bins1=(1:(k1+1))/k1;
% bins1=bins1-(1/k1)/2;
% bins2=(1:(k2+1))/k2;
% bins2=bins2-(1/k2)/2;
bins1=(0:k1)/k1;
bins2=(0:k2)/k2;
%%
figure
set(gcf,'Color','w')
w=6.83;
h=1/3*w;
set(gcf,'Units','inches');
set(gcf,'Position',[11.5 3 w h]);

subplot(1,2,1)
[c1,~]=hist(freqs(1,:),bins1);
c1=c1/sum(c1);
[c2,~]=hist(col(freqs(2:end,:)),bins2);
c2=c2/sum(c2);
% plot(bins2,c2,'-o','LineWidth',lw,'Color',twocols(2,:));
plot(.5*(1+bins2),hygepdf(0:k2,N,m,k2),'-o','LineWidth',lw,'Color',twocols(2,:));
hold on;
% plot(bins1,c1,'-o','LineWidth',lw,'Color',twocols(1,:))
plot(.5*(1+bins1),hygepdf(0:k1,N,m,k1),'-o','LineWidth',lw,'Color',twocols(1,:))
box off
xlabel('Opinion v_i(1)','FontName',fontname,'FontSize',textfontsz)
ylabel('Probability','FontName',fontname,'FontSize',textfontsz)
leg=legend({num2str(k2),num2str(k1)});
legend boxoff
v=get(leg,'Position');
set(leg,'Position',[.2 v(2:4)]) 
set(gca,'FontName',fontname,'FontSize',labfontsz)
set(gca,'xlim',[.45 1])

subplot(1,2,2)
M=max(max(ranks));
width=.4;
[c1,~]=hist(ranks(1,:),1:M);
[c2,~]=hist(ranks(2,:),1:M);
hold on
bar((1:M)+width/2,c2/sum(c2),width,'FaceColor',twocols(2,:),'EdgeColor','k')
bar((1:M)-width/2,c1/sum(c1),width,'FaceColor',twocols(1,:),'EdgeColor','k')
set(gca,'xlim',[0 M+1])
set(gca,'xtick',1:M,'xticklabel',1:M)
legend({num2str(k2),num2str(k1)})
legend boxoff
box off
set(gca,'FontName',fontname,'FontSize',labfontsz)
xlabel('Rank','FontName',fontname,'FontSize',textfontsz)
ylabel('Frequency','FontName',fontname,'FontSize',textfontsz)

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename='/Users/eleanorbrush/Desktop/hypergeometric.pdf';
% print(filename,'-dpdf','-r300')

print('/Users/eleanorbrush/Desktop/test.pdf','-dpdf','-r300')
%%
N=50;
radius=.2;
b=1;
T=3;
strategy=randi([10,19],N,1);
positions=unifrnd(0,1,N,2);
d=squareform(pdist(positions));
receiver=1;
[~,o]=sort(d(receiver,:));
d=d(o,o);
allreceivers=d(receiver,:)<=radius;
M=zeros(N);
for ind=1:N
    [~, order]=sort(d(ind,:));
    neighbors=order(2:strategy(ind)+1);
    M(ind,neighbors)=1/strategy(ind);
end
% M(1:N+1:end)=-1; %sets diagonal equal to -1
% M=M+transpose(M);
M2=M;
props=sum(M2(:,allreceivers==1),2);
M(1:N+1:end)=-sum(M,2);

beta=zeros(N,1);
beta(allreceivers)=b;

A=M-diag(beta);
v=A*beta;
v2=(M+diag(beta))*beta;
v3=(T*eye(N)+T^2/2*A)*beta;
initialconditions=.5*ones(N,1);
% initialconditions=.5*ones(N,1)+normrnd(0,.01,N,1);

[V,Lambda]=eig(A);
[s,o]=sort(real(diag(Lambda)));
V=V(:,o);
Lambda=Lambda(o,o);
w=find(sigfig(diag(Lambda),13)==0,1);
if isempty(w)
    
    part_diag=inv(Lambda)*(expm(T*Lambda)-eye(N));
    part_diag2=T*eye(N)+T^2/2*Lambda;
    part_diag3=part_diag2;
    part_diag3(Lambda<=-1)=1-exp(-T);
    part_diag4=part_diag2;
    k=-2;
    part_diag4(1:N+1:N*sum(allreceivers))=1/k*(exp(T*k)-1);
%     part_diag6=1.5*T*eye(N)+T^2/2*Lambda;
    part_diag6=part_diag2;
    part_diag6(1:N+1:N*sum(allreceivers))=(-1/2*(exp(-2*T)-1));
%     part_diag7=(-2*T*(exp(-2*T)-1)+T)*eye(N)+T^2/2*Lambda;
%     part_diag7=2*(-1/2*(exp(-2*T)-1))*(1+T^2/2*2)*eye(N)+T^2/2*Lambda;
%     part_diag7=(-1/2*(exp(-2*T)-1-(T-T^2)))*eye(N)+T^2/2*Lambda;
    l=Lambda(1,1);
    part_diag7=(1/l*(exp(l*T)-1)+1/2*(-T-T^2/2*(l)))*eye(N)+1/2*(T*eye(N)+T^2/2*Lambda);
    part=V*part_diag*inv(V)*beta;
    part2=V*part_diag2*inv(V)*beta;
    part3=V*part_diag3*inv(V)*beta;
    part4=V*part_diag4*inv(V)*beta;
    part5=V*part_diag2*inv(V)*beta;
    part6=V*part_diag6*inv(V)*beta;
    part7=V*part_diag7*inv(V)*beta;
    particular=expm(T*A)*inv(A)*(eye(N)-expm(-T*A))*beta; %#ok<*MINV>
%     keep=[];
%     Lambda_red=Lambda(keep,keep);
%     V_red=V(:,keep);
%     unkeep=setdiff(1:N,keep);
%     Lambda_left=Lambda(unkeep,unkeep);
%     V_left=V(:,unkeep);
%     particular_reduced=V_red*inv(Lambda_red)*(expm(T*Lambda_red)-eye(length(keep)))*transpose(V_red)*beta;
%     part3=V_red*(T*eye(length(keep))+T^2/2*Lambda_red)*transpose(V_red)*beta;
%     part4=V*inv(Lambda)*(expm(T*Lambda)-eye(N))*transpose(V)*beta;
%     part5=V_left*inv(Lambda_left)*(expm(T*Lambda_left)-eye(length(unkeep)))*transpose(V_left)*beta;
%     part6=V_left*(T*eye(length(unkeep))+T^2/2*Lambda_left)*transpose(V_left)*beta;
else
    fun = @(x)expm(-x*A)*beta;
    integrated = integral(fun,0,T,'ArrayValued',true);
    particular=expm(T*A)*integrated;
end

homogeneous=expm(T*A)*initialconditions;
instantaneous=homogeneous+particular;

check=((1-T/2)*eye(N)+T/2*M2-T/2*diag(beta))*beta;

q=instantaneous;
% [col((corr(particular(setdiff(1:N,find(allreceivers==1))),abs(V(setdiff(1:N,find(allreceivers==1)),:))))); sum(allreceivers)]
% plot(abs(V(setdiff(1:N,find(allreceivers==1)),end)),q(setdiff(1:N,find(allreceivers==1))),'o')
% plot(sum(M2(setdiff(1:N,find(allreceivers==1)),allreceivers==1),2),q(setdiff(1:N,find(allreceivers==1))),'o')
% plot(props(setdiff(1:N,find(allreceivers==1))),particular(setdiff(1:N,find(allreceivers==1))),'o')
% plot(props(setdiff(1:N,find(allreceivers==1))),abs(V(setdiff(1:N,find(allreceivers==1)),end)),'o')
% plot(check(~allreceivers),particular(~allreceivers),'o')

cla
% plot(particular,real(part7),'o');hold on;
% plot(get(gca,'xlim'),(2)*get(gca,'xlim'))

% plot(particular(~allreceivers),props(~allreceivers),'or')
% hold on
% plot(particular(~allreceivers),part5(~allreceivers),'o','MarkerSize',8)
% plot(particular(~allreceivers),part7(~allreceivers),'og')
% s=2*(-1/2*(exp(-2*T)-1))*(1+T^2/2*2)-T;
% plot(get(gca,'xlim'),(s)*get(gca,'xlim'))

plot(particular(),props(),'or')
hold on
plot(particular(),part4(),'ok','MarkerSize',8)
plot(particular(),part5(),'o','MarkerSize',8)
plot(particular(),part7(),'og')
s=2*(-1/2*(exp(-2*T)-1))*(1+T^2/2*2)-T;
plot(get(gca,'xlim'),(s)*get(gca,'xlim'))
plot(get(gca,'xlim'),get(gca,'xlim'),'k')

clc
[corr(particular,props,'type','spearman') corr(particular,real(part2),'type','spearman') corr(particular,real(part3),'type','spearman') corr(particular,real(part4),'type','spearman') corr(particular,real(part5),'type','spearman')    corr(particular,real(part7),'type','spearman')]

%%
p=.1;
N=200;
k=17;
m=floor(p*N);

combin=zeros(min(k,m)+1,1);
binom=zeros(min(k,m)+1,1);
hypo=zeros(min(k,m)+1,1);

for j=0:min(k,m)
    combin(j+1)=nchoosek(k,j)*nchoosek(N-k,m-j)/nchoosek(N,m);
%     binom(j+1)=nchoosek(k,j)*power(p,j)*power(1-p,k-j);
    hypo(j+1)=hygepdf(j,N,m,k);
end

%%
N=50;
strategy=47*ones(N,1);
strategy(1)=45;
m=9;
k=15;
toadd=zeros(m+1,1);
for l=0:m
   toadd(l+1)=hygepdf(l,N-1,m,strategy(1))*power(1-hygecdf(l,N-1,m,strategy(1)),k-1)*power(1-hygecdf(l/strategy(1)*strategy(2),N-1,m,strategy(2)),N-k);
end
probless=sum(toadd);
%%
N=20;
resident=17;
invader=13;
m=9;
k=1;

toadd=zeros(m+1,1);
for l=0:m
   toadd(l+1)=hygepdf(l,N-1,m,invader)*power(1-hygecdf(l,N-1,m,invader),k-1)*power(1-hygecdf(floor(l/invader*resident),N-1,m,resident),N-k);
end
probleast=sum(toadd);
toadd2=zeros(m+1,1);
for l=0:m
   toadd2(l+1)=hygepdf(l,N-1,m,resident)*power(1-hygecdf(l,N-1,m,resident),N-k-1)*power(1-hygecdf(floor(l/resident*invader),N-1,m,invader),k);
end
probleast2=sum(toadd2);
probeaten=zeros(1,N);
probeaten(1:k)=probleast;
probeaten((k+1):end)=(1-probleast*k)/(N-k);

toadd3=zeros(m+1,1);
for l=0:m
%                toadd(l+1)=hygepdf(l,N-1,m,invader)*power(hygecdf(l-1,N-1,m,invader),k-1)*power(hygecdf(l/invader*resident,N-1,m,resident),N-k);
   toadd3(l+1)=hygepdf(l,N-1,m,invader)*power(hygecdf(l-1,N-1,m,invader),k-1)*power(hygecdf((l-1)/invader*resident,N-1,m,resident),N-k); %l-1 because (l-1)/invader*resident isn't going to be an integer anyway
end
probmost=sum(toadd3);
probgettoeat=zeros(1,N);
probgettoeat(1:k)=probmost;
probgettoeat((k+1):end)=(1-probmost*k)/(N-k);
%%
N=20;
resident=1;
invader=1;
m=9;
k=1;
toadd=zeros(m+1,1);
for l=0:m
   toadd(l+1)=hygepdf(l,N-1,m,invader)*power(1-hygecdf(l,N-1,m,invader),k-1)*power(1-hygecdf(floor(l/invader*resident),N-1,m,resident),N-k);
   toadd(l+1)=toadd(l+1)+hygepdf(l,N-1,m,invader)*sum(power(hygepdf(l,N-1,m,invader),1:(k-1))./(2:k))*power(1-hygecdf(floor(l/invader*resident),N-1,m,resident),N-k);
end
probleast=sum(toadd);
probeaten=zeros(1,N);
probeaten(1:k)=probleast;
probeaten((k+1):end)=(1-probleast*k)/(N-k);

toadd2=zeros(m+1,1);
for l=0:m
%                toadd(l+1)=hygepdf(l,N-1,m,invader)*power(hygecdf(l-1,N-1,m,invader),k-1)*power(hygecdf(l/invader*resident,N-1,m,resident),N-k);
   toadd2(l+1)=hygepdf(l,N-1,m,invader)*power(hygecdf(l-1,N-1,m,invader),k-1)*power(hygecdf(floor(l/invader*resident),N-1,m,resident),N-k);
end
probmost=sum(toadd2);
probgettoeat=zeros(1,N);
probgettoeat(1:k)=probmost;
probgettoeat((k+1):end)=(1-probmost*k)/(N-k);

