function M=makenet_normed(strategy)
strategy=reshape(strategy,1,length(strategy));
N=size(strategy);
N=N(2);

positions=zeros(N,2);
positions(:,1)=cos((1:N)*2*pi/N);
positions(:,2)=sin((1:N)*2*pi/N);
d=squareform(pdist(positions));

M=zeros(N);
    for ind=1:N
        deg=strategy(1,ind);
        if strategy(1,ind)>1 && strategy(1,ind)<N
        [~, order]=sort(d(ind,:));
        neighbors=order(2:strategy(1,ind)+1);
        M(ind,neighbors)=1/deg;
        end
        if strategy(1,ind)==1
            if ind==N
            M(ind,1)=1/deg;
            else M(ind,ind+1)=1/deg;
            end
        end
        if strategy(1,ind)==N
            M(ind,:)=1/deg;
            M(ind,ind)=0;
        end
    end
   M(1:N+1:end)=-1;
end
    