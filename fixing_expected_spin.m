% strategy=randi(19,20,1);
strategy=[19,5*ones(1,19)];
positions=unifrnd(0,1,N,2);
d=squareform(pdist(positions));
k=3;
beta=[ones(k,1);zeros(N-k,1)];

M=zeros(N);
for ind=1:N
    [~, order]=sort(d(ind,:));
    neighbors=order(2:strategy(ind)+1);
    M(ind,neighbors)=1/strategy(ind);
end
M(1:N+1:end)=-1;
A=M-diag(beta);

[vecs,vals]=eig(A);
V=vecs;
vals=diag(vals);
nonzero=abs(vals)>0.00001;
nonzerovals=vals(nonzero);

%if there is a zero eigenvalue set the eigenvector to all ones
%if there are two eigenvectors for eigenvalue 0 I specify what they should
%be here
indices=1:N;
where=indices(~nonzero);
vhat=1/sqrt(N)*ones(N,1);
if sum(~nonzero)==1
    V(:,where(1))=ones(N,1);
%     V(:,where(1))=vhat;
end
if sum(~nonzero)>1
    L=sum(~nonzero);
    for i=1:(L-1)
        vproj=V(:,where(i+1))-(V(:,where(i+1))'*vhat)*vhat;
        V(:,where(i+1))=vproj;
    end
end

% mu=pinv(V)*beta;s

% c=zeros(N,1);
% c(where(1))=1;

initialconditions=.5*ones(N,1);
c=pinv(V)*initialconditions;

l=exp(nonzerovals*T);
Lhom=zeros(1,N);
Lhom(nonzero)=l;
Lhom(~nonzero)=1;
Lhom=diag(Lhom);
Phi=V*Lhom;
Uhom=Phi;

l=1./nonzerovals.*(exp(nonzerovals*T)-1);
Linhom=zeros(1,N);
Linhom(nonzero)=l;
Linhom(~nonzero)=T;
Linhom=diag(Linhom);
% Uinhom=V*Linhom/V;
Uinhom=V*Linhom*pinv(V);

instantaneous=Uinhom*beta+Uhom*c;

% l=1./power(nonzerovals,2).*(exp(nonzerovals*T)-1)-T./nonzerovals;
% L=zeros(1,N);
% L(nonzero)=l;
% L(~nonzero)=T^2/2;
% L=diag(L);
% U=V*L/V;
% integrated=U*beta;

% q=[instantaneous,integrated];

q=instantaneous;


[V2,vals2]=eig(V);
vals2=diag(vals2);
indices=1:N;
nonzero=indices(sigfig(vals2,3)~=0);
zero=indices(sigfig(vals2,3)==0);
k=1;
R=V2(:,nonzero(k));
% [V*pinv(V)*R-R]
inv(V);
q=V*pinv(V)*beta-beta;
% max(abs(q))

%%
nowfun = @(x) A*x;
[t,y] = ode45(@sqrt,[0 1],.5);
