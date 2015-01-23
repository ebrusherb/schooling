strategy=(N-1)*ones(N,1);
positions=unifrnd(0,1,N,2);
d=squareform(pdist(positions));

M=zeros(N);
for ind=1:N
    [~, order]=sort(d(ind,:));
    neighbors=order(2:strategy(ind)+1);
    M(ind,neighbors)=1/strategy(ind);
end
M(1:N+1:end)=-1;

k=0;

q=zeros(N-1,N);
q(1,1:2)=[1/sqrt(2) -1/sqrt(2)];
for j=2:(N-1)
   v=rand(1,N);
   v=v-sum(v)/N;
   for k=1:(j-1)
      v=v-(v*transpose(q(k,:)))*q(k,:); 
   end
   v=v/norm(v,2);
   q(j,:)=v;
end


A=M-diag(beta);
P=-q*A*transpose(q);
Pinv=inv(P);

[vecs,vals]=eig(P);
vals=diag(vals);
invvals=1./vals;
U=vecs*diag(invvals)*transpose(vecs);

sigma1=q*transpose(A)*ones(N,1);
sigma2=q*beta;

zcorrs=zeros(N-1,N-1);
for i=1:(N-1)
    for j=i:(N-1)
%         piece1=sum(vecs(i,:).*invvals'.*vecs(j,:));
%         piece1=sum(vecs(i,:).*reshape(invvals,1,[]).*vecs(j,:));
%         piece2=sum(sum((sigma1'*Pinv(:,i))*(sigma1'*Pinv(:,j))'))+sum(sum((sigma1'*Pinv(:,j))*(sigma1'*Pinv(:,i))'));
%         piece3=sum(sum((sigma2'*Pinv(:,i))*(sigma2'*Pinv(:,j))'))+sum(sum((sigma2'*Pinv(:,j))*(sigma2'*Pinv(:,i))'));
%         piece2=(reshape(sigma1,1,[])*Pinv(:,i))*(reshape(sigma1,1,[])*col(Pinv(j,:)));
%         piece3=(reshape(sigma2,1,[])*Pinv(:,i))*(reshape(sigma2,1,[])*col(Pinv(j,:)));
        Pinv_deriv=kron(Pinv(:,i),Pinv(j,:));
        piece1=U(i,j);
        piece2=reshape(sigma1,1,[])*Pinv_deriv*col(sigma1);
        piece3=reshape(sigma2,1,[])*Pinv_deriv*col(sigma2);
        zcorrs(i,j)=1/2*piece1+1/4*piece2+1/4*piece3;
%         zcorrs(j,i)=1/2*piece1+1/4*piece2+1/4*piece3;
%         zcorrs(i,j)=piece1+piece2+piece3;
%         zcorrs(j,i)=piece1+piece2+piece3;
    end
end

% corrs=transpose(q)*zcorrs*q;

cov=transpose(q)*zcorrs*q;
todivide=repmat(diag(cov),1,N).*repmat(diag(cov)',N,1);
todivide=power(todivide,.5);
corrs=cov./todivide;