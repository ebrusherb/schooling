% matlabpool local 2

numsigs_permove=1;
nummoves=1000;
N=20;
b=1;

strats=2:2:(N-1);
L=length(strats);
 
T=1;
radius=.5;

fitnesseaten=zeros(L,L);
rhoeaten=zeros(L,L);

fitnessgettoeat=zeros(L,L);
rhogettoeat=zeros(L,L);

featen=zeros(L,L,N-1);
geaten=zeros(L,L,N-1);

fgettoeat=zeros(L,L,N-1);
ggettoeat=zeros(L,L,N-1);

tic
parfor ind=1:(L*L*(N-1))
    [u,v,k]=ind2sub([L,L,N-1],ind);
    squareind=sub2ind([L,L],u,v);
    resident=strats(u);
    invader=strats(v);

    strategy=resident*ones(1,N);
    strategy(1:k)=invader;
    [meanscore, scorevar, probeaten, probgettoeat]=signalingevents(strategy,numsigs_permove,nummoves,radius,b,T);

    perfeaten=1-probeaten;
    featen(ind)=mean(perfeaten(1:k));
    geaten(ind)=mean(perfeaten((k+1):N));

    perfgettoeat=probgettoeat;
    fgettoeat(ind)=mean(perfgettoeat(1:k));
    ggettoeat(ind)=mean(perfgettoeat((k+1):N));
         
%     if k==N-1
%         fitnesseaten(squareind)=featen(k)/geaten(k);
%         fitnessgettoeat(squareind)=fgettoeat(k)/ggettoeat(k);
%     end
%     ratio=geaten./featen;
%     tosum=ones(1,N-1);
%     for k=1:(N-1)
%         tosum(k)=prod(ratio(1:k));
%     end
%     rhoeaten(squareind)=1/(1+sum(tosum));
% 
%     ratio=ggettoeat./fgettoeat;
%     tosum=ones(1,N-1);
%     for k=1:(N-1)
%         tosum(k)=prod(ratio(1:k));
%     end
%     rhogettoeat(squareind)=1/(1+sum(tosum));
%     end
end
toc

tic
for ind=1:L*L
    [u,v]=ind2sub([L,L],ind);
    where=zeros(N-1,1);
    for k=1:N-1
        where(k)=sub2ind([L,L,N-1],u,v,k);
    end
    featen_now=featen(where);
    geaten_now=geaten(where);
    
    fgettoeat_now=fgettoeat(where);
    ggettoeat_now=ggettoeat(where);
    
    fitnesseaten(ind)=featen_now(1)/geaten_now(1);
    fitnessgettoeat(ind)=fgettoeat_now(1)/ggettoeat_now(1);
    
    ratio=featen_now/geaten_now;
    tosum=ones(1,N-1);
    for k=1:(N-1)
        tosum(k)=prod(ratio(1:k));
    end
    rhoeaten(ind)=1/(1+sum(tosum));
    
    ratio=fgettoeat_now/ggettoeat_now;
    tosum=ones(1,N-1);
    for k=1:(N-1)
        tosum(k)=prod(ratio(1:k));
    end
    rhogettoeat(ind)=1/(1+sum(tosum));
end   
toc
matlabpool close
%%
% matlabpool local 2

its=20;

N=40;

stratvals=[5 10 15 20];
bvals=[0 0.01 1 2 5 7 10];
radvals=[.01 .1 .5 1];

Ns=length(stratvals);
Nb=length(bvals);
Nr=length(radvals);

corrmat=zeros(Nb,its);
H2mat=zeros(its,1);

for itnow=1:its
    positions=unifrnd(0,1,N,2);
    d=squareform(pdist(positions));

    corrvec=zeros(Nb,1);

    strat=stratvals(3);
    strategy=strat*ones(N,1);

    radius=radvals(3);

    M=zeros(N);
    for ind=1:N
        [~, order]=sort(d(ind,:));
        neighbors=order(2:strategy(ind)+1);
        M(ind,neighbors)=1/strategy(ind);
    end
    M(1:N+1:end)=-1;
    H2mat(itnow,1)=H2norm(M);
    tic    
    parfor j=1:Nb 
        b=bvals(j);
        [distbins,avgcorr,corrlength] = correlationlength_mat(M,d,b,radius);
        corrvec(j)=corrlength;
    end
    toc
corrmat(:,itnow)=corrvec;

end   
%%
maxvals=max(corrmat);
[r,c]=find(corrmat==repmat(maxvals,Nb,1));
minb=accumarray(c,r,[],@min);
minb=bvals(minb);
% matlabpool close

