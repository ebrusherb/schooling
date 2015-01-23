
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

noise=eye(N);
noise=q*noise*transpose(q);

L=-M;

Lbar=q*L*transpose(q);

sigma = lyap(Lbar,-noise);

beta=zeros(N,1);

A=M-diag(beta);
P=-q*A*transpose(q);
Pinv=inv(P);

[vecs,vals]=eig(P);
vals=diag(vals);
invvals=1./vals;

sigma1=q*transpose(A)*ones(N,1);
sigma2=q*beta;

zcorrs=zeros(N-1,N-1);
for i=1:(N-1)
    for j=i:(N-1)
%         piece1=sum(vecs(i,:).*invvals'.*vecs(j,:));
%         piece1=sum(vecs(i,:).*reshape(invvals,1,[]).*vecs(j,:));
        piece1=U(i,j);
        piece2=sum(sum((sigma1'*Pinv(:,i))*(sigma1'*Pinv(:,j))'))+sum(sum((sigma1'*Pinv(:,j))*(sigma1'*Pinv(:,i))'));
        piece3=sum(sum((sigma2'*Pinv(:,i))*(sigma2'*Pinv(:,j))'))+sum(sum((sigma2'*Pinv(:,j))*(sigma2'*Pinv(:,i))'));
        zcorrs(i,j)=piece1+piece2+piece3;
        zcorrs(j,i)=piece1+piece2+piece3;
    end
end

