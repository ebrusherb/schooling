N=20;
strategy=[5,5*ones(1,N-1)];
positions=unifrnd(0,1,N,2);
d=squareform(pdist(positions));

M=zeros(N);
for ind=1:N
    [~, order]=sort(d(ind,:));
    neighbors=order(2:strategy(ind)+1);
    M(ind,neighbors)=1/strategy(ind);
end
M(1:N+1:end)=-1;

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
%%
k=3;
beta=[ones(k,1);zeros(N-k,1)];

A=M-diag(beta);

L=M;

noise=eye(N);

Lbar=q*L*transpose(q);
noise=q*noise*transpose(q);
sigma = lyap(Lbar,noise);
h = sqrt(trace(sigma));

y0=.5*ones(N,1);
R=q*beta;

fun1 = @(x)expm((T-x)*Lbar)*R*transpose(R);
integrated1 = integral(fun1,0,T,'ArrayValued',true);
fun2 = @(x)R*transpose(R)*expm((T-x)*Lbar);
integrated2 = integral(fun2,0,T,'ArrayValued',true);

sigma_forced = lyap(Lbar,noise+integrated1+integrated2);
h_forced = sqrt(trace(sigma_forced));
[h h_forced]

M1=max(max(abs(sigma)));
M2=max(max(abs(sigma_forced)));
Mtot=max(M1,M2);
subplot(1,2,1);
imagesc(sigma);
caxis manual
caxis([-Mtot Mtot])
colorbar;
subplot(1,2,2);
imagesc(sigma_forced);
caxis manual
caxis([-Mtot Mtot])
colorbar

%%
% fun1 = @(x)expm((T-x)*A)*beta*transpose(beta);
% integrated1 = integral(fun1,0,T,'ArrayValued',true);
% fun2 = @(x)beta*transpose(beta)*expm((T-x)*A);
% integrated2 = integral(fun2,0,T,'ArrayValued',true);

% sigma_forced = lyap(A,eye(N)+integrated1+integrated2);

sigma_forced = lyap(A,eye(N)-inv(A)*beta*transpose(beta)-beta*transpose(beta)*transpose(inv(A)));
delta = lyap(transpose(A),eye(N));

%%
s1=lyap(A,eye(N));
V=sigma_forced-inv(A)*beta*transpose(beta)*transpose(inv(A));
s2=lyap(transpose(A),eye(N));
Q=transpose(A)*inv(V)+inv(V)*A;
s3=lyap(transpose(A),-Q);
% Q2=(transpose(A)*A+transpose(A)*transpose(A)+A*A+transpose(A)*A);
% Q2=A+transpose(A);
% Q2=-2*diag([16*ones(k,1);4*ones(N-k,1)]);
% Q2=-diag(diag(Q));
% s4=lyap(transpose(A),-Q2);
s4 = -inv(A+transpose(A));

% s1=lyap(A,eye(N)+integrated1+integrated2);
% s2=lyap(A,eye(N)-inv(A)*beta*transpose(beta)-beta*transpose(beta)*inv(A));
s=cell(4,1);
s{1}=s1;
s{2}=s2;
s{3}=s3;
s{4}=s4;
subplot(1,5,1)
imagesc(s{1})
colorbar
for i=2:4
    subplot(1,5,i)
    imagesc(inv(s{i}))
    colorbar
end
subplot(1,5,5)
imagesc(V);
colorbar
%%
T=1000;
fun1 = @(x)expm((T-x)*A)*beta*transpose(beta);
integrated1 = integral(fun1,0,T,'ArrayValued',true);

% check=inv(A)*(expm(T*A)-eye(N))*beta*transpose(beta);
check=inv(A)*(-eye(N))*beta*transpose(beta);

subplot(1,2,1)
imagesc(integrated1)
colorbar
subplot(1,2,2)
imagesc(check)
colorbar

%%
sigma_f=lyap(A,eye(N));
sigma_f_transpose=lyap(transpose(A),eye(N));
Q=transpose(A)*inv(sigma_f)+inv(sigma_f)*A;
Q2=A+transpose(A);
M=A+transpose(A);
Q3=(transpose(A)+A)*(M)+(M)*(transpose(A)+A);
% Q3=transpose(A)*M+M*A;
C=lyap(transpose(A),eye(N));
C1=lyap(transpose(A),-Q);
C2=lyap(transpose(A),-Q2);
C3=inv(-Q2);
dynamics=lyap((M),eye(N));

subplot(2,3,1)
imagesc(sigma_f);colorbar
subplot(2,3,2)
imagesc(C);colorbar
subplot(2,3,3)
imagesc(C1);colorbar
subplot(2,3,4)
imagesc(C2);colorbar
subplot(2,3,5)
% imagesc(C3);colorbar
imagesc(dynamics);colorbar
subplot(2,3,6)
imagesc(sigma_f_transpose);colorbar


%%
N=50;
strategy=[15,15*ones(1,N-1)];
positions=unifrnd(0,1,N,2);
d=squareform(pdist(positions));

k=0;
b=0.1;
beta=b*[ones(k,1);zeros(N-k,1)];

L=zeros(N);
M=zeros(N);
for ind=1:N
    [~, order]=sort(d(ind,:));
    neighbors=order(2:strategy(ind)+1);
    M(ind,neighbors)=1/strategy(ind);
    L(ind,neighbors)=1/strategy(ind);
end
M(1:N+1:end)=-1;

Nmat=.5*(L+transpose(L));
Nbar=Nmat-diag(Nmat*ones(N,1));
B=diag(beta);
P=Nbar-B;
A=M-B;
% 
% sigma=lyap(A,eye(N));
% 
% D=lyap(P,eye(N));
% 
% figure
% subplot(1,4,1)
% imagesc(A)
% colorbar
% subplot(1,4,2)
% imagesc(D)
% colorbar
% subplot(1,4,3)
% imagesc(sigma)
% colorbar
% subplot(1,4,4)
% imagesc(-inv(P))
% colorbar
% 
% Q=A*inv(P)+inv(P)*transpose(A);
% 
% todivide=repmat(diag(-inv(P)),1,N).*repmat(reshape(diag(-inv(P)),1,[]),N,1);
% todivide=power(todivide,.5);
% corrs=-inv(P)./todivide;
% 
% figure
% plot(d,corrs,'ok')
% set(gca,'ylim',[0 1])
% 
% [vecs,vals]=eig(P);

[W,Lambda]=eig(Nbar);
w=find(sigfig(diag(Lambda),10)==0);
Wtilde=W(:,[(1:w-1) (w+1):N]);
Lambdatilde=Lambda([(1:w-1) (w+1):N],[(1:w-1) (w+1):N]);
Ntilde=Wtilde*Lambdatilde*transpose(Wtilde);
Ntildeinv=Wtilde*inv(Lambdatilde)*transpose(Wtilde);

todivide=repmat(diag(-Ntildeinv),1,N).*repmat(reshape(diag(-Ntildeinv),1,[]),N,1);
todivide=power(todivide,.5);
corrs=-Ntildeinv./todivide;

dcol=col(d);
corrscol=col(corrs);
[~,o]=sort(dcol);
plot(dcol(o),corrscol(o),'ok')
hold on
plot(get(gca,'xlim'),zeros(2,1),'r')
% hold off
set(gca,'ylim',[-1 1])

%%
Ntest=[-1 1;1 -1];
[W,Lambda]=eig(Ntest);
w=find(sigfig(diag(Lambda),10)~=0);
Wtilde=W(:,w);
Lambdatilde=Lambda(w,w);
Ntilde=Wtilde*Lambdatilde*transpose(Wtilde);
Ntildeinv=Wtilde*inv(Lambdatilde)*transpose(Wtilde);
