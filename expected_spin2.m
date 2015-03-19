function q = expected_spin2(M,T,beta)
%solves the differential equation dq/dt=M*q+beta*(1-q)  -->
%dq/dt=A*q+beta where A=M-diag(beta)

% N=length(strategy);
% M=makenet_normed(strategy);
A=M-diag(beta);
N=size(M,1);

if nargin==2
    beta=1*ones(N,1);
end

initialconditions=.5*ones(N,1);

[~,vals]=eig(A);
w=find(sigfig(diag(vals),13)==0,1);
if isempty(w)
    particular=expm(T*A)*inv(A)*[eye(N)-expm(-T*A)]*beta;
else
    fun = @(x)expm(-x*A)*beta;
    integrated = integral(fun,0,T,'ArrayValued',true);
    particular=expm(T*A)*integrated;
end

homogeneous=expm(T*A)*initialconditions;
instantaneous=homogeneous+particular;

q=instantaneous;
end

