function [corrs]=paircorrelations(M,beta)

N=size(M,2);

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

% Mstar=q*M*transpose(q);
% gamma=-transpose(M)*ones(N,1)-beta;
% gammastar=q*gamma;
% Mstarinv=inv(Mstar);
% 
% [vecs,vals]=eig(Mstar);
% vals=diag(vals);
% invvals=1./vals;
% 
% zcorrs=zeros(N-1,N-1);
% for i=1:(N-1)
%     for j=i:(N-1)
%         soc=-sum(vecs(i,:).*invvals'.*vecs(j,:));
%         sig=sum(sum((gammastar'*Mstarinv(:,i))*(gammastar'*Mstarinv(:,j))'))+sum(sum((gammastar'*Mstarinv(:,j))*(gammastar'*Mstarinv(:,i))'));
% %         sig=0;
%         zcorrs(i,j)=soc+sig;
%         zcorrs(j,i)=soc+sig;
%     end
% end

negM=-M;
Mstar=q*negM*transpose(q);
Mstarinv=inv(Mstar);

[vecs,vals]=eig(Mstar);
vals=diag(vals);
invvals=1./vals;

sigma1=q*transpose(negM)*ones(N,1);
sigma2=q*beta;

zcorrs=zeros(N-1,N-1);
for i=1:(N-1)
    for j=i:(N-1)
        piece1=sum(vecs(i,:).*invvals'.*vecs(j,:));
        piece2=sum(sum((sigma1'*Mstarinv(:,i))*(sigma1'*Mstarinv(:,j))'))+sum(sum((sigma1'*Mstarinv(:,j))*(sigma1'*Mstarinv(:,i))'));
        piece3=sum(sum((sigma2'*Mstarinv(:,i))*(sigma2'*Mstarinv(:,j))'))+sum(sum((sigma2'*Mstarinv(:,j))*(sigma2'*Mstarinv(:,i))'));
        zcorrs(i,j)=piece1+piece2+piece3;
        zcorrs(j,i)=piece1+piece2+piece3;
    end
end

corrs=transpose(q)*zcorrs*q;
end