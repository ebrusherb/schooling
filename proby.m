function p = proby(y,Pf,Pftilde)
N=length(y);
eN=[zeros(N-1,1);1];
% p=1/sqrt((2*pi)^(N-1)*det(inv(-Pf))*Pftilde(N,N))*exp(-1/2*transpose(y+(V-1)*sqrt(N)*eN)*(Pftilde)*(y+(V-1)*sqrt(N)*eN));
% p=1/sqrt((2*pi)^(N-1)*det(inv(-Pf))*Pftilde(N,N))*exp(-1/2*transpose(y)*Pftilde*y)*exp(-1/2*(Pftilde(N,N)*N*(V-1)^2+2*(V-1)*sqrt(N)*sum(Pftilde(N,:)*y)));
% p=1/sqrt((2*pi)^(N-1)*det(inv(-Pf))*Pftilde(N,N))*exp(-1/2*transpose(y)*(Pftilde-1/Pftilde(N,N)*Pftilde(:,end)*Pftilde(end,:))*y)*exp(-1/2*Pftilde(N,N)*N*(V-1+sum(Pftilde(end,:)*y)/Pftilde(N,N)/sqrt(N))^2);
% p=sqrt(N)/sqrt((2*pi)^N*det(inv(-Pf)))*exp(-1/2*(transpose(y)*(Pftilde-1/Pftilde(N,N)*Pftilde(:,end)*Pftilde(end,:))*y));
p=1/sqrt((2*pi)^(N-1)*det(inv(-Pf))*Pftilde(N,N))*exp(-1/2*(transpose(y)*(Pftilde-1/Pftilde(N,N)*Pftilde(:,end)*Pftilde(end,:))*y));
end