function [analytical numerical toprint M beta]=check_solver(T)

N=evalin('base','N');
strategy=evalin('base','strategy');
radius=evalin('base','radius');
b=evalin('base','b');

positions=unifrnd(0,1,N,2);
d=squareform(pdist(positions));

M=zeros(N);
for ind=1:N
    [~, order]=sort(d(ind,:));
    neighbors=order(2:strategy(ind)+1);
    M(ind,neighbors)=1/strategy(ind);
end
M(1:N+1:end)=-1; %sets diagonal equal to -1
receiver=randsample(N,1,'true');
beta=zeros(N,1);
allreceivers=d(receiver,:)<=radius;
beta(allreceivers)=b;
analytical=expected_spin(M,T,beta);

if max(abs(imag(analytical)))>0.001
    toprint='problem';
else toprint='nope';
end

k=100;
deltat=T/k;
q=ones(N,1);
for i=1:k
    q=(M*q+beta)*deltat+q;
end

numerical=q;
    
end