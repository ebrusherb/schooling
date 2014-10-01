Nb1=10;
Nb2=10;
l=1;
keeptrackofinfoD=zeros(its,Nb1*Nb2);
keeptrackofinfoR=zeros(its,Nb1*Nb2);
keeptrackofinfoDelta=zeros(its,Nb1*Nb2);
keeptrackofinfoH=zeros(its,Nb1*Nb2);
keeptrackofinfoPi=zeros(its,Nb1*Nb2);
infomatD=zeros(Nb1,Nb2);
infomatR=zeros(Nb1,Nb2);
infomatDelta=zeros(Nb1,Nb2);
infomatH=zeros(Nb1,Nb2);
infomatPi=zeros(Nb1,Nb2);

powerfuns={'D','R','Delta','H','Pi'};
numpowerfuns=length(powerfuns);

skewvalsD=zeros(1,its);
skewvalsR=zeros(1,its);
skewvalsDelta=zeros(1,its);
skewvalsH=zeros(1,its);
skewvalsPi=zeros(1,its);
finalthresh=zeros(N,its);

keepprobmat=zeros(its,N*N);
keeptimemat=zeros(its,N*N);
keepsigmat=zeros(its,N*N);

tic
for c=1:its

    opt_its=10;

    fighting_abilities=zeros(1,N);
    
    switch dist
        case 'unif'
            fighting_abilities=20*rand(1,N);
        case 'norm'
            fighting_abilities=normrnd(0,10,1,N);
        case 'lognorm'
            fighting_abilities=lognrnd(0,1,1,N);
    end
    
    fighting_abilities=sort(fighting_abilities,'descend');

    dmat=zeros(N); %turn difference in fighting abilities into probabilitise of winning
    extendeddomvals=[1-domvals(end:-1:2) domvals];
    v=1:(Nd+Nd-1);
    deps=.05;

    for i=1:N
        for j=[1:(i-1),(i+1):N]
            diff=fighting_abilities(i)-fighting_abilities(j);
            d=exp(diff)/(exp(diff)+1);
            d=floor(round(d/deps))*deps;
            dmat(i,j)=v(abs(d-extendeddomvals)<=deps/2);
        end
    end

    %size(twomat)=[Nl,Nd,Nt,Nt,2];
    
    perf=zeros(2,Nd,Nt,Nt); %individual, leak, dominance, thresholds

    perf(1,2:Nd,:,:)=c1*(1-twomat(l,2:Nd,:,:,1))+c2*(twomat(l,2:Nd,:,:,2))+c3*(1-twomat(l,2:Nd,:,:,1));
    perf(2,2:Nd,:,:)=c1*(1-twomat(l,2:Nd,:,:,1))+c2*(twomat(l,2:Nd,:,:,2))+c3*(twomat(l,2:Nd,:,:,1));

    perf(1,1,:,:)=c2*(twomat(l,1,:,:,2))+c3*(1-twomat(l,1,:,:,1));
    perf(2,1,:,:)=c2*(twomat(l,1,:,:,2))+c3*(twomat(l,1,:,:,1));

    Tvals=zeros(N,opt_its);
    Tvals(:,1)=2*ones(N,1);
    % Tvals(:,1)=randi([1 Nt],N,1);
    perfvals=zeros(N,opt_its);

    for i=1:N
        opp_thresh=Tvals([1:(i-1),(i+1):N],1);
        ds=dmat(i,[1:(i-1),(i+1):N]);
        perfsum=0;
                for q=1:(N-1)
                    if ds(q)<=Nd-1
                        perfsum=perfsum+perf(2,Nd-ds(q)+1,opp_thresh(q),Tvals(i,1)); 
                    else 
                        perfsum=perfsum+perf(1,ds(q)-Nd+1,Tvals(i,1),opp_thresh(q));
                    end
                end  
        perfvals(i,1)=perfsum;
    end

    count=1;
    while count<=opt_its
        for i=1:N
            opp_thresh=Tvals([1:(i-1),(i+1):N],count);
            ds=dmat(i,[1:(i-1),(i+1):N]);
            perftest=zeros(1,Nt);
            for j=1:Nt
                perfsum=0;
                for q=1:(N-1)
                    if ds(q)<=Nd-1
                        perfsum=perfsum+perf(2,Nd-ds(q)+1,opp_thresh(q),j);
                    else 
                        perfsum=perfsum+perf(1,ds(q)-Nd+1,j,opp_thresh(q));
                    end
                end 
                perftest(j)=perfsum;
            end
            [p,n]=min(perftest);
            Tvals(i,count+1)=n;
            perfvals(i,count+1)=p;
        end
        if sum(Tvals(:,count+1)==Tvals(:,count))==N
            maxit=count+2;
            count=opt_its+1;
            Tvals(:,maxit:end)=[];
            perfvals(:,maxit:end)=[];
        end
        count=count+1;
    end
    
    finalthresh(:,c)=threshvals(Tvals(:,end));
    if c==its
        threshtohist=threshvals(Tvals(:,end));
    end

    probmat=zeros(N,N);
    timemat=zeros(N,N);

    %size(twomat)=[Nl,Nd,Nt,Nt,2];

    for i=1:N
        for j=[1:(i-1),(i+1):N]
            if dmat(i,j)<=Nd-1
                probmat(i,j)=1-twomat(l,Nd-dmat(i,j)+1,Tvals(j,end),Tvals(i,end),1);
                timemat(i,j)=twomat(l,Nd-dmat(i,j)+1,Tvals(j,end),Tvals(i,end),2);
            else
                probmat(i,j)=twomat(l,dmat(i,j)-Nd+1,Tvals(i,end),Tvals(j,end),1);
                timemat(i,j)=twomat(l,dmat(i,j)-Nd+1,Tvals(i,end),Tvals(j,end),2);
            end
        end
    end
    
    binsigmat=zeros(N,N);
    
    for i=1:N
        for j=(i+1):N
            draw=rand;
            if draw <+probmat(i,j)
                binsigmat(i,j)=1;
            else binsigmat(j,i)=1;
            end
        end
    end
    
    sigmat=binsigmat.*(max(max(timemat))+1-timemat);
    
    todivide=sum(sigmat,2);
    todivide(todivide==0)=1;
    entmat=sigmat./repmat(todivide,1,N);
    logmat=entmat;
    logmat(entmat~=0)=-log(entmat(entmat~=0));
    entmat=entmat.*logmat;

%         powervec=sum(probmat,2);
    powervecD=sum(binsigmat,2);
    powervecR=sum(sigmat,2);
    powervecDelta=(powervecD).*(powervecR);
    powervecH=sum(entmat,2);
    powervecH=real(powervecH);
    powervecPi=(powervecR).*(powervecH);
    
    abilitiesdiff=ceil(max(fighting_abilities))-floor(min(fighting_abilities));
    abilitiesbins=floor(min(fighting_abilities)):abilitiesdiff/Nb1:ceil(max(fighting_abilities));
    [~,abilitiesindices]=histc(fighting_abilities,abilitiesbins);
    
    powervec=powervecD;
    powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
    powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
    [~,powerindices]=histc(powervec,powerbins);
    
    v=zeros(Nb1*Nb2,1);
    for i=1:N
        ind=sub2ind([Nb1,Nb2],abilitiesindices(i),powerindices(i));
        v(ind)=v(ind)+1;
    end
    keeptrackofinfoD(c,:)=v;
    
    skewvalsD(c)=skewness(powervec);
    
    powervec=powervecR;
    powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
    powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
    [~,powerindices]=histc(powervec,powerbins);
    
    v=zeros(Nb1*Nb2,1);
    for i=1:N
        ind=sub2ind([Nb1,Nb2],abilitiesindices(i),powerindices(i));
        v(ind)=v(ind)+1;
    end
    keeptrackofinfoR(c,:)=v;
    
    skewvalsR(c)=skewness(powervec);
    
    powervec=powervecDelta;
    powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
    powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
    [~,powerindices]=histc(powervec,powerbins);
    
    v=zeros(Nb1*Nb2,1);
    for i=1:N
        ind=sub2ind([Nb1,Nb2],abilitiesindices(i),powerindices(i));
        v(ind)=v(ind)+1;
    end
    keeptrackofinfoDelta(c,:)=v;
    
    skewvalsDelta(c)=skewness(powervec);
    
    powervec=powervecH;
    powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
    powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
    [~,powerindices]=histc(powervec,powerbins);
    
    v=zeros(Nb1*Nb2,1);
    for i=1:N
        ind=sub2ind([Nb1,Nb2],abilitiesindices(i),powerindices(i));
        v(ind)=v(ind)+1;
    end
    keeptrackofinfoH(c,:)=v;
    
    skewvalsH(c)=skewness(powervec);
    
    powervec=powervecPi;
    powerdiff=ceil(max(powervec))-floor(min(powervec))+.1;
    powerbins=(floor(min(powervec))-.05):(powerdiff/Nb2):(ceil(max(powervec))+.05);
    [~,powerindices]=histc(powervec,powerbins);
    
    v=zeros(Nb1*Nb2,1);
    for i=1:N
        ind=sub2ind([Nb1,Nb2],abilitiesindices(i),powerindices(i));
        v(ind)=v(ind)+1;
    end
    keeptrackofinfoPi(c,:)=v;
    
    skewvalsPi(c)=skewness(powervec);
 
        keepprobmat(c,:)=reshape(probmat,1,[]);
        keeptimemat(c,:)=reshape(timemat,1,[]);
        keepsigmat(c,:)=reshape(binsigmat,1,[]);
    
end
toc
keeptrackofinfoD=sum(keeptrackofinfoD,1);
keeptrackofinfoR=sum(keeptrackofinfoR,1);
keeptrackofinfoDelta=sum(keeptrackofinfoDelta,1);
keeptrackofinfoH=sum(keeptrackofinfoH,1);
keeptrackofinfoPi=sum(keeptrackofinfoPi,1);

for i=1:Nb1
    for j=1:Nb2
        ind=sub2ind([Nb1,Nb2],i,j);
        infomatD(i,j)=keeptrackofinfoD(ind);
    end
end
for i=1:Nb1
    for j=1:Nb2
        ind=sub2ind([Nb1,Nb2],i,j);
        infomatR(i,j)=keeptrackofinfoR(ind);
    end
end
for i=1:Nb1
    for j=1:Nb2
        ind=sub2ind([Nb1,Nb2],i,j);
        infomatDelta(i,j)=keeptrackofinfoDelta(ind);
    end
end
for i=1:Nb1
    for j=1:Nb2
        ind=sub2ind([Nb1,Nb2],i,j);
        infomatH(i,j)=keeptrackofinfoH(ind);
    end
end
for i=1:Nb1
    for j=1:Nb2
        ind=sub2ind([Nb1,Nb2],i,j);
        infomatPi(i,j)=keeptrackofinfoPi(ind);
    end
end

avgthresh=mean(finalthresh,2);

probmat=reshape(keepprobmat(its,:),N,N);
timemat=reshape(keeptimemat(its,:),N,N);
binsigmat=reshape(keepsigmat(its,:),N,N);

groupmutinfoD=0;

for i=1:Nb1
    px=sum(infomatD(i,:))/(N*its);
    for j=1:Nb2
        py=sum(infomatD(:,j))/(N*its);
        pxy=infomatD(i,j)/(N*its);
        if pxy~=0
            groupmutinfoD=groupmutinfoD+pxy*log(pxy/px/py);
        end
    end
end

groupmutinfoR=0;

for i=1:Nb1
    px=sum(infomatR(i,:))/(N*its);
    for j=1:Nb2
        py=sum(infomatR(:,j))/(N*its);
        pxy=infomatR(i,j)/(N*its);
        if pxy~=0
            groupmutinfoR=groupmutinfoR+pxy*log(pxy/px/py);
        end
    end
end

groupmutinfoDelta=0;

for i=1:Nb1
    px=sum(infomatDelta(i,:))/(N*its);
    for j=1:Nb2
        py=sum(infomatDelta(:,j))/(N*its);
        pxy=infomatDelta(i,j)/(N*its);
        if pxy~=0
            groupmutinfoDelta=groupmutinfoDelta+pxy*log(pxy/px/py);
        end
    end
end
groupmutinfoH=0;

for i=1:Nb1
    px=sum(infomatH(i,:))/(N*its);
    for j=1:Nb2
        py=sum(infomatH(:,j))/(N*its);
        pxy=infomatH(i,j)/(N*its);
        if pxy~=0
            groupmutinfoH=groupmutinfoH+pxy*log(pxy/px/py);
        end
    end
end
groupmutinfoPi=0;

for i=1:Nb1
    px=sum(infomatPi(i,:))/(N*its);
    for j=1:Nb2
        py=sum(infomatPi(:,j))/(N*its);
        pxy=infomatPi(i,j)/(N*its);
        if pxy~=0
            groupmutinfoPi=groupmutinfoPi+pxy*log(pxy/px/py);
        end
    end
end

meanskewnessD=mean(skewvalsD);
meanskewnessR=mean(skewvalsR);
meanskewnessDelta=mean(skewvalsDelta);
meanskewnessH=mean(skewvalsH);
meanskewnessPi=mean(skewvalsPi);

% toreturn={avgthresh,probmat,timemat,sigmat,groupmutinfo,meanskewness};