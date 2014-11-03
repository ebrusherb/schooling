function q = expected_spin(M,T,beta)
%solves the differential equation dq/dt=M*q+beta*(1-q)  -->
%dq/dt=A*q+beta where A=M-diag(beta)

% N=length(strategy);
% M=makenet_normed(strategy);
A=M-diag(beta);
N=size(M,1);

if nargin==2
    beta=1*ones(N,1);
end

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
end

