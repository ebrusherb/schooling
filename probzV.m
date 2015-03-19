function p = probzV(z,V,Pf)
N=length(z);
Qfull=evalin('base','Qfull');
y=Qfull*z;
eN=[zeros(N-1,1);1];
% p=1/sqrt((2*pi)^N*det(inv(-Pf)))*exp(-1/2*transpose(z-(1-V)*ones(N,1))*(-Pf)*(z-(1-V)*ones(N,1)));
% p=exp(-1/2*transpose(v-ones(N,1))*(-Pf)*(v-ones(N,1)));
p=1/sqrt((2*pi)^N*det(inv(-Pf)))*exp(-1/2*transpose(transpose(Qfull)*y-(1-V)*sqrt(N)*transpose(Qfull)*eN)*(-Pf)*(transpose(Qfull)*y-(1-V)*sqrt(N)*transpose(Qfull)*eN));
end