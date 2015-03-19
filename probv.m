function p = probv(v,Pf)
N=length(v);
p=1/sqrt((2*pi)^N*det(inv(-Pf)))*exp(-1/2*transpose(v-ones(N,1))*(-Pf)*(v-ones(N,1)));
% p=exp(-1/2*transpose(v-ones(N,1))*(-Pf)*(v-ones(N,1)));
end