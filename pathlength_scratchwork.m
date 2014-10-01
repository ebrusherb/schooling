N=10;
strategy=4*ones(N,1);
positions=unifrnd(0,1,N,2);
    d=squareform(pdist(positions));

    M=zeros(N);
    for ind=1:N
        [~, order]=sort(d(ind,:));
        neighbors=order(2:strategy(ind)+1);
        M(ind,neighbors)=1;
    end
    g=sparse(M);
    
    paths=graphallshortestpaths(g);
    upperpaths=triu(paths,1);
    lowerpaths=triu(transpose(paths),1);
    pathsvec=[upperpaths(~~upperpaths), lowerpaths(~~lowerpaths)]
    pathsvec=min(pathsvec')';
    
    [distbins,avgcorr,corrlength] = correlationlength(M,d,b,radius);
    
    