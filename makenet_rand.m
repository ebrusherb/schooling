function M=makenet_rand(strategy)
N=max(size(strategy));

positions=unifrnd(0,1,N,2);
d=squareform(pdist(positions));

M=zeros(N);
    for ind=1:N
%         if strategy(ind)>1 && strategy(ind)<N
        [~, order]=sort(d(ind,:));
        neighbors=order(2:strategy(1,ind)+1);
        M(ind,neighbors)=1;
%         end
%         if strategy(ind)==1
%             if ind==N
%             M(ind,1)=1;
%             else M(ind,ind+1)=1;
%             end
%         end
%         if strategy(ind)==N
%             M(ind,:)=1;
%             M(ind,ind)=0;
%         end
    end
end
    