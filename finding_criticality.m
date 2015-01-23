
T=1;
teststrategy=2*randi([1 floor((N-1)/2)],N,1);
positions=unifrnd(0,1,N,2);
d=squareform(pdist(positions));
[~,receiver]=min(teststrategy(:,end));
radius=1.5;
bvals=[0 100 500 100000 1000000];
Nb=length(bvals);
storeddistbins=cell(Nb,1);
storedavgcorr=cell(Nb,1);
corrlengths=zeros(Nb,1);


M=zeros(N);
for ind=1:N
    [~, order]=sort(d(ind,:));
    neighbors=order(2:teststrategy(ind)+1);
    M(ind,neighbors)=1/teststrategy(ind);
end
M(1:N+1:end)=-1;
    
    for j=1:Nb 
        b=bvals(j);
%             [distvec,corrvec,corrlength]=correlationlength_mat(M,d,b,radius);
        [distvec,corrvec,corrlength] = correlationlength_mat_single(M,d,b,radius,receiver);
        
        storeddistbins{j}=distvec;
        storedavgcorr{j}=corrvec;
        corrlengths(j)=corrlength;
    end
    
    %%
figure 
hold on

for i=1:Nb
plot(storeddistbins{i},storedavgcorr{i})
end
i=Nb;plot(storeddistbins{i},storedavgcorr{i},'-or')
    
%%
N=20;
T=1;
% strategy=2*randi([1 floor((N-1)/2)],N,1);
strategy=8*ones(N,1);
strategy(end)=2;
positions=unifrnd(0,1,N,2);
d=squareform(pdist(positions));
receiver=1;
% [~,o]=sort(d(receiver,:),'ascend');
% d=d(o,o);

b=1;
radvals=0:.05:1.6;
% radvals=sort([radvals .7],'ascend');
Nr=length(radvals);
storeddistvecs=cell(Nr,1);
storedcorrvecs=cell(Nr,1);
corrlengths=zeros(Nr,1);
numrecs=zeros(Nr,1);
scores=zeros(Nr,N);

M=zeros(N);
for ind=1:N
    [~, order]=sort(d(ind,:));
    neighbors=order(2:strategy(ind)+1);
    M(ind,neighbors)=1/strategy(ind);
end
M(1:N+1:end)=-1;
    
    for j=1:Nr
        radius=radvals(j);
%         [distvec,corrvec,corrlength] = correlationlength_mat(M,d,b,radius);
        [distvec,corrvec,corrlength,~] = correlationlength_mat_single_fromrec(M,d,b,radius,receiver);
%         [distvec,corrvec,corrlength]=correlationlength_mat_single(M,d,b,radius,receiver);
        distvec(isnan(corrvec))=[];
        corrvec(isnan(corrvec))=[];
        storeddistvecs{j}=distvec;
        storedcorrvecs{j}=corrvec;
        corrlengths(j)=corrlength;
        numrecs(j)=sum(d(receiver,:)<=radius);
        beta=zeros(N,1);
        beta(d(receiver,:)<=radius)=b;
        scores(j,:)=expected_spin(M,T,beta);
    end
    
radlabs=cell(Nr,1);
for i=1:Nr
    radlabs{i}=num2str(radvals(i));
end
%%
close all
% mycols=transpose(fliplr(transpose(cbrewer('seq','PuBuGn',Nr))));
mycols=colormap(parula(Nr));
figure 

subplot(2,3,1:3)
hold on
plot(get(gca,'xlim'),zeros(2,1),'k')
p=zeros(1,Nr);
for i=1:Nr
    p(i)=plot(storeddistvecs{i},storedcorrvecs{i},'-o','Color',mycols(i,:),'LineWidth',2);
end
legend(p,radlabs)

    
subplot(2,3,4)
hold on
plot(radvals,corrlengths,'-o')
plot(radvals,numrecs/N,'-or')
plot(radvals,radvals,'-k')
plot(radvals,.5*ones(Nr,1),'-k')

subplot(2,3,5)
plot(numrecs/N,corrlengths,'-o')
set(gcf,'Position',[40 378 1300 700])

subplot(2,3,6)
hold on
for i=1:Nr
    plot(storedcorrvecs{i},scores(i,:)/scores(i,receiver),'o','Color',mycols(i,:),'LineWidth',2)
end
