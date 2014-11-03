N=20;
radius=1.4;
T=.1;
nummoves=1000;
numsigs_permove=1;
b=0;

strategy=[17,16*ones(1,N-1)];
numsigs_tot=numsigs_permove*nummoves;
N=max(size(strategy));
indvec=1:N;

positions=unifrnd(0,1,N,2);
d=squareform(pdist(positions));

M=zeros(N);
for ind=1:N
    [~, order]=sort(d(ind,:));
    neighbors=order(2:strategy(ind)+1);
    M(ind,neighbors)=1/strategy(ind);
end
M(1:N+1:end)=-1; %sets diagonal equal to -1
receivers=randsample(N,numsigs_permove,'true');
j=1;
beta=zeros(N,1);
receiver=receivers(j);
allreceivers=d(receiver,:)<=radius;
beta(allreceivers)=b;
v=real(expected_spin(M,T,beta));
v2=real(expected_spin(M,T,beta));

perm=randsample(20,20,'false');
Mperm=M(perm,perm);
betaperm=beta(perm);
vperm=real(expected_spin(Mperm,T,betaperm));
% [v(perm) vperm v2]

%%
A=M-diag(beta);
N=size(M,1);

[vecs,vals]=eig(A);
V=vecs;
vals=diag(vals);
nonzero=abs(vals)>0.00001;
nonzerovals=vals(nonzero);

Aperm=Mperm-diag(betaperm);

[vecsperm,valsperm]=eig(Aperm);
Vperm=vecsperm;
valsperm=diag(valsperm);
nonzeroperm=abs(valsperm)>0.00001;
nonzerovalsperm=valsperm(nonzeroperm);

indices=1:N;
vhat=1/sqrt(N)*ones(N,1);
where=indices(~nonzero);
if sum(~nonzero)==1
    V(:,where(1))=ones(N,1);
end
if sum(~nonzero)>1
    L=sum(~nonzero);
    for i=1:(L-1)
        vproj=V(:,where(i+1))-(V(:,where(i+1))'*vhat)*vhat;
        V(:,where(i+1))=vproj;
    end
end

whereperm=indices(~nonzeroperm);
if sum(~nonzeroperm)==1
    Vperm(:,whereperm(1))=ones(N,1);
end
if sum(~nonzeroperm)>1
    L=sum(~nonzeroperm);
    for i=1:(L-1)
        vprojperm=Vperm(:,whereperm(i+1))-(Vperm(:,whereperm(i+1))'*vhat)*vhat;
        Vperm(:,whereperm(i+1))=vprojperm;
    end
end
    
d1=diag(V);
d2=diag(Vperm);
% [d1(perm) d2]

%%
radius=1.5;
T=1;
nummoves=1000;
numsigs_permove=1;
b=1;
perm=randsample(20,20,'false');

strategy=[16*ones(1,N)];
strategy(7)=17;
% strategy(7)=19;
numsigs_tot=numsigs_permove*nummoves;
N=max(size(strategy));
indvec=1:N;

scores = zeros(N, numsigs_permove, nummoves);
keepreceivers = zeros(N, numsigs_permove, nummoves);
keepvals = zeros(N,nummoves);
keepvecs = zeros(N,N,nummoves);

% parfor i=1:nummoves
for i = 1:nummoves
    positions=unifrnd(0,1,N,2);
    d=squareform(pdist(positions));

    M=zeros(N);
    for ind=1:N
        [~, order]=sort(d(ind,:));
        neighbors=order(2:strategy(ind)+1);
        M(ind,neighbors)=1/strategy(ind);
    end
    M(1:N+1:end)=-1; %sets diagonal equal to -1
    
    receivers=randsample(N,numsigs_permove,'true');
    j=1;
        beta=zeros(N,1);
        receiver=receivers(j);
        allreceivers=d(receiver,:)<=radius;
        beta(allreceivers)=b;
        v=real(expected_spin(M,T,beta));
       
        scores(:,j,i)=v;
        keepreceivers(:,j,i)=(allreceivers);
    
    A=M-diag(beta);
    N=size(M,1);

    [vecs,vals]=eig(A);
    V=vecs;
    vals=diag(vals);
    [~,o]=sort(real(vals));
    temp=V(:,o);
%     s=sign(temp(1,:));
%     temp=temp*diag(s);
    keepvecs(:,:,i)=temp';
    keepvals(:,i)=vals(o);
end
scores=reshape(scores,N,[]);
stdevscores=std(scores,[],2);

reshaped=reshape(keepvecs,400,1000);
summed=reshape(sum(keepvecs,1),N,[]);
resstd=std(reshaped,[],2);
sumstd=std(summed,[],2);
matsd=reshape(resstd,20,[]);
matsd=transpose(matsd);

subplot(1,3,1)
plot(sum(matsd,2),stdevscores,'o')

subplot(1,3,2)
plot(sumstd,stdevscores,'o')

subplot(1,3,3)
plot(sumstd)

%%
strats=1:19;
L=length(strats);
lambdavals=zeros(L,1);
H2vals=zeros(1,L);
corrlengthvals=zeros(1,L);
discvals=zeros(1,L);

numsigs_permove=1;
nummoves=1000;
b=1;
radius=1;
T=1;


for i=1:L
    strategy=strats(i)*ones(N,1);
    [meanlambda, meanH2, meancorrlength, disconnectedcount]=groupprops(strategy,numsigs_permove,nummoves,radius,b,T);
    lambdavals(i)=meanlambda;
    H2vals(i)=meanH2;
    corrlengthvals(i)=meancorrlength;
    discvals(i)=disconnectedcount;
end

%%
strat=18;
strategy=strat*ones(N,1);

positions=unifrnd(0,1,N,2);
d=squareform(pdist(positions));

M=zeros(N);
for ind=1:N
    [~, order]=sort(d(ind,:));
    neighbors=order(2:strategy(ind)+1);
    M(ind,neighbors)=1/strategy(ind);
end
M(1:N+1:end)=-1; %sets diagonal equal to -1
receivers=randsample(N,numsigs_permove,'true');
j=1;
beta=zeros(N,1);
receiver=receivers(j);
allreceivers=d(receiver,:)<=radius;
beta(allreceivers)=b;

A=M-diag(beta);
N=size(M,1);

[vecs,vals]=eig(A);
V=vecs;
vals=diag(vals);

x=sqrt((N-2)^2/((N-2)^2+(N-1)));
y=-1/(N-2)*x;
lam=-(N-1)/(N-2);

testmat=y*ones(N,N);
testmat(1:N+1:N^2)=x;
check=A*testmat-testmat*diag(lam);
sumcheck=sum(check,1);
weird=sigfig(sumcheck,10)~=0;

Mbin=M;
Mbin(M==-1)=0;
Mbin(~~Mbin)=1;
g=sparse(Mbin);
[~,Ctotal]=graphconncomp(g,'Directed','false');
[paths] = graphallshortestpaths(g);

gplotdc(M,positions)
rainbow=colormap(jet(max(Ctotal)));
hold on
for i=1:max(Ctotal)
    plot(positions(Ctotal==i,1),positions(Ctotal==i,2),'o','LineWidth',5,'Color',rainbow(i,:))
end
plot(positions(allreceivers,1),positions(allreceivers,2),'ok','MarkerSize',10,'LineWidth',3)
plot(positions(weird,1),positions(weird,2),'or','MarkerSize',15,'LineWidth',3)
colorbar
hold off

q = expected_spin(M,T,beta)
