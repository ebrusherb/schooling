function [h,s] = networkprops(M,scaled,beta)
L=-lap(M);
N=size(M,1);

s=size(L);
n=s(1);
if nargin==1
    scaled='additive';
end
if strcmp(scaled,'additive')==1
noise=diag(sum(transpose(M)));
else
noise=eye(n);
end
L(1:N+1:end)=diag(L)+beta;
sigma = lyap(L,noise);
h = (trace(sigma))/N;
s = sum(sum(sigma))/N^2;
end