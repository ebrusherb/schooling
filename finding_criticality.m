% [distvec,corrvec,corrlength,crosscorrs] = correlationlength_mat_single_v3(M,d,b,radius,i)
N=50;
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
%         [distvec,corrvec,corrlength,~] = correlationlength_mat_single_v3(M,d,b,radius,receiver);
        [distvec,corrvec,corrlength]=correlationlength_mat(M,d,b,radius);
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
    