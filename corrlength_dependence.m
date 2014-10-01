%%

its=20;

N=40;

stratvals=[10 15 20];
bvals=[0 0.01 .1:.1:1 1.5];
radvals=[.01 .1 .5 1];

Ns=length(stratvals);
Nb=length(bvals);
Nr=length(radvals);

%%
%find correlation length, correlation between correlation function at b and
%at b=1000, and H2 norm for networks with various strategies and radii of
%signal
corrlengthmat=zeros(its,Ns,Nr);
similaritymat=zeros(its,Ns,Nr,Nb,2);
minbmat=zeros(its,Ns,Nr);
H2mat=zeros(its,Ns);
for itnow=1:its
    positions=unifrnd(0,1,N,2);
    d=squareform(pdist(positions));

    for i=1:Ns
        strat=stratvals(i);
        strategy=strat*ones(N,1);

        M=zeros(N);
        for ind=1:N
            [~, order]=sort(d(ind,:));
            neighbors=order(2:strategy(ind)+1);
            M(ind,neighbors)=1/strategy(ind);
        end
        M(1:N+1:end)=-1;
        H2mat(itnow,i)=H2norm(M);
        
        for k=1:Nr
            radius=radvals(k);
            similarityvec=zeros(Nb,2);
            [~,endavgcorr,corrlength] = correlationlength_mat(M,d,1000,radius);
            stop1=find(isnan(endavgcorr),1,'first')-1;
            corrlengthmat(itnow,i,k)=corrlength;
            [~,beginavgcorr,~] = correlationlength_mat(M,d,00,radius);
            stop2=find(isnan(beginavgcorr),1,'first')-1;
            stop=min(stop1,stop2);
            parfor j=1:Nb 
                b=bvals(j);
                [~,avgcorr,corrlength] = correlationlength_mat(M,d,b,radius);
                similarityvec(j,:)=[corr(avgcorr(1:stop)',beginavgcorr(1:stop)'), corr(avgcorr(1:stop)',endavgcorr(1:stop)')]; %#ok<PFBNS>
            end
            similaritymat(itnow,i,k,:,:)=similarityvec;
            v=similarityvec(:,1)-similarityvec(:,2);
            boundary=sum(v>0);
            minbmat(itnow,i,k)=interp1([v(boundary) v(boundary+1)],[bvals(boundary) bvals(boundary+1)],0);
        end
    end
end 
%%
%plot correlation between correlation function as a function of strategy
%for various matrices
itnow=itnow+1;plot(bvals,reshape(similaritymat(itnow,:,k,:),Ns,[]),'-o');set(gca,'xlim',[0 2])
%%
%plot correlation functions for various b and 2 different strategies to
%convince myself that as strategy increases correlation length increases
positions=unifrnd(0,1,N,2);
d=squareform(pdist(positions));

bvals2=[0 0.01 .1 .5 1 5 10:10:50];
Nb2=length(bvals2);
storeddistbins=cell(Nb2,2);
storedavgcorr=cell(Nb2,2);
ref=cell(1,2);
corrlengths=zeros(1,2);
for i=[1 3]
    
    strat=stratvals(i);
    strategy=strat*ones(N,1);
    
    M=zeros(N);
    for ind=1:N
        [~, order]=sort(d(ind,:));
        neighbors=order(2:strategy(ind)+1);
        M(ind,neighbors)=1/strategy(ind);
    end
    M(1:N+1:end)=-1;
    
    radius=.5;
    distbinsvecs=cell(Nb2,1);
    avgcorrvecs=cell(Nb2,1);
    [distbins,avgcorr,corrlength] = correlationlength_mat(M,d,10000,radius);
    ref{ceil(i/2)}=avgcorr;
    corrlengths(ceil(i/2))=corrlength;
    parfor j=1:Nb2 
        b=bvals(j);
        [distbins,avgcorr,corrlength] = correlationlength_mat(M,d,b,radius);
        distbinsvecs{j}=distbins;
        avgcorrvecs{j}=avgcorr;
    end
    
    for j=1:Nb2
        storeddistbins{j,ceil(i/2)}=distbinsvecs{j};
        storedavgcorr{j,ceil(i/2)}=avgcorrvecs{j};
    end
    
end
%%
c=colormap(jet(Nb2));
figure
for i=1:2
    subplot(1,2,i)
    hold on
    for j=1:Nb2
        plot(storeddistbins{j,i},storedavgcorr{j,i},'Color',c(j,:),'LineWidth',2)
    end
    plot(storeddistbins{Nb2,i},ref{i},'k','LineWidth',2)
    plot([0 1.5],zeros(2,1))
end

compcorrmat=zeros(2,Nb2);
for i=1:2
    for j=1:Nb2
        compcorrmat(i,j)=corr(ref{i}(1:50)',storedavgcorr{j,i}(1:50)');
    end
end
        
    
%%
strategy1=10*ones(N,1);

strategy2=20*ones(N,1);

b=1000;

diff=0;
count=0;
while diff<=0 && count<=100
    positions=unifrnd(0,1,N,2);
    d=squareform(pdist(positions));
    M1=zeros(N);
    for ind=1:N
        [~, order]=sort(d(ind,:));
        neighbors=order(2:strategy1(ind)+1);
        M1(ind,neighbors)=1/strategy1(ind);
    end
    M1(1:N+1:end)=-1;
    M2=zeros(N);
    for ind=1:N
        [~, order]=sort(d(ind,:));
        neighbors=order(2:strategy2(ind)+1);
        M2(ind,neighbors)=1/strategy2(ind);o
    end
    M2(1:N+1:end)=-1;

    [distbins,avgcorr1,corrlength1] = correlationlength_mat(M1,d,b,radius);  
    [distbins,avgcorr2,corrlength2] = correlationlength_mat(M2,d,b,radius);
    diff=corrlength1-corrlength2;
    count=count+1;
end

%%
%plot correlation function for a single receiver as a function of strategy
%of network
N=40;
radius=.5;
b=10;

positions=unifrnd(0,1,N,2);
d=squareform(pdist(positions));

receiver=1;
beta=zeros(N,1);
allreceivers=d(receiver,:)<=radius;
beta(allreceivers)=b;

strategy1=10*ones(N,1);

M1=zeros(N);
for ind=1:N
    [~, order]=sort(d(ind,:));
    neighbors=order(2:strategy1(ind)+1);
    M1(ind,neighbors)=1/strategy1(ind);
end
M1(1:N+1:end)=-1;

corrs1=paircorrelations(M1,beta);

uppercorr1=triu(corrs1,1);
upperdist1=triu(d,1);

corrvec1=[uppercorr1(~~uppercorr1); diag(corrs1)];
distvec1=[upperdist1(~~upperdist1); diag(d)];
[s,o]=sort(distvec1);

distvec1=distvec1(o);
corrvec1=corrvec1(o);

strategy2=30*ones(N,1);

M2=zeros(N);
for ind=1:N
    [~, order]=sort(d(ind,:));
    neighbors=order(2:strategy2(ind)+1);
    M2(ind,neighbors)=1/strategy2(ind);
end
M2(1:N+1:end)=-1;

receiver=1;
beta=zeros(N,1);
allreceivers=d(receiver,:)<=radius;
beta(allreceivers)=b;

corrs2=paircorrelations(M2,beta);

uppercorr2=triu(corrs2,1);
upperdist2=triu(d,1);

corrvec2=[uppercorr2(~~uppercorr2); diag(corrs2)];
distvec2=[upperdist2(~~upperdist2); diag(d)];
[s,o]=sort(distvec2);

distvec2=distvec2(o);
corrvec2=corrvec2(o);

%%
its=10;
N=100;

stratvals=[10 20];
bvals=[0 .5 1 1.5 5 10 100];
radvals=[.1 1.5];

Ns=length(stratvals);
Nb=length(bvals);
Nr=length(radvals);

receiver=1;
distmat=zeros(its,N);
corrfunmat=zeros(its,Ns,Nr,Nb,N);
corrlengthmat=zeros(its,Ns,Nr,Nb);
similaritymat=zeros(its,Ns,Nr,Nb,2);
minbmat=zeros(its,Ns,Nr);
H2mat=zeros(its,Ns);
for itnow=1:its
    positions=unifrnd(0,1,N,2);
    d=squareform(pdist(positions));
    distmat(itnow,:)=d(receiver,:);
    
    for i=1:Ns
        strat=stratvals(i);
        strategy=strat*ones(N,1);

        M=zeros(N);
        for ind=1:N
            [~, order]=sort(d(ind,:));
            neighbors=order(2:strategy(ind)+1);
            M(ind,neighbors)=1/strategy(ind);
        end
        M(1:N+1:end)=-1;
        H2mat(itnow,i)=H2norm(M);
        
        for k=1:Nr
            radius=radvals(k);
            similarityvec=zeros(Nb,2);
            [~,endcorrvec,corrlength] = correlationlength_mat_single_v2(M,d,1000,radius,receiver);
            [~,begincorrvec,~] = correlationlength_mat_single_v2(M,d,00,radius,receiver);
            valid=~isnan(begincorrvec);
            storecorrfun=zeros(Nb,N);
            storecorrlength=zeros(Nb,1);
            for j=1:Nb 
                b=bvals(j);
                [~,corrvec,corrlength] = correlationlength_mat_single_v2(M,d,b,radius,receiver);
                storecorrfun(j,:)=corrvec;
                storecorrlength(j)=corrlength;
                similarityvec(j,:)=[corr(corrvec(valid)',begincorrvec(valid)'), corr(corrvec(valid)',endcorrvec(valid)')]; %#ok<PFBNS>
            end
            corrfunmat(itnow,i,k,:,:)=storecorrfun;
            corrlengthmat(itnow,i,k,:)=storecorrlength;
%             similaritymat(itnow,i,k,:,:)=similarityvec;
%             v=similarityvec(:,1)-similarityvec(:,2);
%             boundary=sum(v>0);
%             if boundary<Nb
%                 minbmat(itnow,i,k)=interp1([v(boundary) v(boundary+1)],[bvals(boundary) bvals(boundary+1)],0);
%             else minbmat(itnow,i,k)=bvals(end);
%             end
        end
    end
end 

%%
it=1;
distbins=0:.02:1.5;

i=1;
k=2;

figure
c=colormap(jet(Nb));
% hold on
for i=1:2
    subplot_tight(1,2,i,[.05 .05])
    hold on
% [~,o]=sort(distmat(it,:));
for j=1:Nb
    plot(sort(distmat(it,:)),reshape(corrfunmat(it,i,k,j,:),1,[]),'Color',c(j,:),'LineWidth',2);
end
plot(distmat(it,:),zeros(N,1),'k')
% plot(radvals(k)*ones(2,1),get(gca,'ylim'),'k');
end
set(gcf,'Position',[600  378 700 600])
%%
it=1;
maxdiff=zeros(Ns,Nr,Nb);

for i=1:Ns
    for k=1:Nr
        for j=1:Nb
            maxdiff(i,k,j)=max(corrfunmat(it,i,k,j,:))-min(corrfunmat(it,i,k,j,:));
        end
    end
end
%%
it=it+1;
divide=find(sort(distmat(it,:))>radvals(k),1,'first');
colors=['r' 'b'];
figure
hold on
for i=1:Ns
    firstdiff=corrfunmat(it,i,k,1,divide-1)-corrfunmat(it,i,k,1,divide);
    plot(bvals,(reshape(corrfunmat(it,i,k,:,divide-1)-corrfunmat(it,i,k,:,divide),1,[]))/firstdiff,colors(i))
end
set(gcf,'Position',[600  378 700 600])

%%
its=20;
N=100;

stratvals=[10 20];
bvals=[0 .5 1 1.5 5 10 100];
radvals=[.1];

Ns=length(stratvals);
Nb=length(bvals);
Nr=length(radvals);

receiver=1;

crosscorrmat=zeros(its,Ns,Nr,Nb,4);
H2mat=zeros(its,Ns);

for itnow=1:its
    positions=unifrnd(0,1,N,2);
    d=squareform(pdist(positions));
    distmat(itnow,:)=d(receiver,:);
    
    for i=1:Ns
        strat=stratvals(i);
        strategy=strat*ones(N,1);

        M=zeros(N);
        for ind=1:N
            [~, order]=sort(d(ind,:));
            neighbors=order(2:strategy(ind)+1);
            M(ind,neighbors)=1/strategy(ind);
        end
        M(1:N+1:end)=-1;
        H2mat(itnow,i)=H2norm(M);
        
        for k=1:Nr
            radius=radvals(k);
            
            storecrosscorr=zeros(Nb,4);
            for j=1:Nb 
                b=bvals(j);
                [~,~,~,crosscorrs] = correlationlength_mat_single_v3(M,d,b,radius,receiver);
                storecrosscorr(j,:)=crosscorrs;
            end
            crosscorrmat(itnow,i,k,:,:)=storecrosscorr;
            
        end
    end
end 

%%
it=it+1;
i=1;
k=1;

figure
c=colormap(jet(Nb));
for i=1:Ns
subplot_tight(1,Ns,i,[.05 .05])
hold on
% for j=1:Nb
%     plot(reshape(crosscorrmat(it,i,k,j,:),1,[]),'-o','Color',c(j,:))
% end
% plot(bvals,reshape(crosscorrmat(it,i,k,:,4),Nb,[]),'-o')
plot(bvals,reshape(crosscorrmat(it,i,k,:,3),Nb,[])./abs(repmat(reshape(crosscorrmat(it,i,k,1,3),1,[]),Nb,1)),'-o')
end

% v1=reshape(crosscorrmat(it,1,k,:,3),Nb,[])./abs(repmat(reshape(crosscorrmat(it,1,k,1,3),1,[]),Nb,1));
% v2=reshape(crosscorrmat(it,2,k,:,3),Nb,[])./abs(repmat(reshape(crosscorrmat(it,2,k,1,3),1,[]),Nb,1));
% plot(v1,v2)

set(gcf,'Position',[600  378 700 600])

%%
v1=reshape(crosscorrmat(:,1,k,end,3),its,[])./abs(reshape(crosscorrmat(:,1,k,1,3),its,1));
v2=reshape(crosscorrmat(:,2,k,end,3),its,[])./abs(reshape(crosscorrmat(:,2,k,1,3),its,1))