function [meanlambda, meanH2, meancorrlength, disconnectedcount]=groupprops_parallel(strategy,numsigs_permove,nummoves,radius,b,T)

N=max(size(strategy));

lambdavals=zeros(1,nummoves);
H2norms=zeros(1,nummoves);
corrlengths=zeros(numsigs_permove,nummoves);

disconnectedcount=0;

parfor pari=1:nummoves
% for pari=1:nummoves
     
    positions=unifrnd(0,1,N,2);
    d=squareform(pdist(positions));

    M=zeros(N);
    for ind=1:N
        [~, order]=sort(d(ind,:));
        neighbors=order(2:(strategy(ind)+1));
        M(ind,neighbors)=1/strategy(ind);
    end
    M(1:N+1:end)=-1; %sets diagonal equal to -1
        
    Mbin=M;
    Mbin(M==-1)=0;
    Mbin(~~Mbin)=1;
    g=sparse(Mbin);
    [~,Ctotal]=graphconncomp(g,'Directed','false');

    L=lap(M);
    [~,vals]=eig(L);
    vals=diag(vals);
    f=find(sigfig(vals,14)==0);

    if max(Ctotal)==1 && length(f)==1           
        vals=vals(sigfig(vals,14)~=0);
        lambda=min(real(vals));
        lambdavals(pari)=lambda;
        h=H2norm(M,'additive');
        H2norms(pari)=h; 
    else 
        disconnectedcount=disconnectedcount+1;
        lambdavals(pari)=NaN;
        H2norms(pari)=NaN;
    end
    
    receivers=randsample(N,numsigs_permove,'true');
    for j=1:numsigs_permove
        receiver=receivers(j);
        
        if max(Ctotal)==1 && length(f)==1
%             [~,~,l,~]=correlationlength_mat_single_v3(M,d,b,radius,receiver);
            [~,~,l] = correlationlength_mat_single(M,d,b,radius,receiver);
            corrlengths(j,pari)=l;
        else 
            corrlengths(j,pari)=NaN;
        end
    end
end

meanlambda=mean(lambdavals(~isnan(lambdavals)));
meanH2=mean(H2norms(~isnan(H2norms)));
meancorrlength=mean(col(corrlengths(~isnan(corrlengths))));
