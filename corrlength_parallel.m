its=10;

N=40;

stratvals=[10 15 20];
% stratvals=[10 20];
bvals=[0 0.01 .1 .5 .75 1:5];
radvals=[.01 .1 .5 1];
% radvals=.5;

Ns=length(stratvals);
Nb=length(bvals);
Nr=length(radvals);


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
            corrlengthmat(itnow,i,k)=corrlength;
            [~,beginavgcorr,~] = correlationlength_mat(M,d,00,radius);
            valid=~isnan(beginavgcorr);
            parfor j=1:Nb 
                b=bvals(j);
                [~,avgcorr,corrlength] = correlationlength_mat(M,d,b,radius);
                similarityvec(j,:)=[corr(avgcorr(valid)',beginavgcorr(valid)'), corr(avgcorr(valid)',endavgcorr(valid)')]; %#ok<PFBNS>
            end
            similaritymat(itnow,i,k,:,:)=similarityvec;
            v=similarityvec(:,1)-similarityvec(:,2);
            boundary=sum(v>0);
            if boundary==Nb
                minbmat(itnow,i,k)=bvals(end);
            else 
                minbmat(itnow,i,k)=interp1([v(boundary) v(boundary+1)],[bvals(boundary) bvals(boundary+1)],0);
            end
        end
    end
end 



N=40;

stratvals=[10 15 20];
radvals=[.01 .1 .5 1];

Ns=length(stratvals);
Nb=length(bvals);
Nr=length(radvals);

positions=unifrnd(0,1,N,2);
d=squareform(pdist(positions));

bvals2=[0 0.01 .1 .5 1 5 10:10:50 1000 10000 50000];
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
    [distbins,avgcorr,corrlength] = correlationlength_mat(M,d,1000,radius);
    ref{ceil(i/2)}=avgcorr;
    corrlengths(ceil(i/2))=corrlength;
    parfor j=1:Nb2 
        b=bvals2(j);
        [distbins,avgcorr,corrlength] = correlationlength_mat(M,d,b,radius);
        distbinsvecs{j}=distbins;
        avgcorrvecs{j}=avgcorr;
    end
    
    for j=1:Nb2
        storeddistbins{j,ceil(i/2)}=distbinsvecs{j};
        storedavgcorr{j,ceil(i/2)}=avgcorrvecs{j};
    end
    
end