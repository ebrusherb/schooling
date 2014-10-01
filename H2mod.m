function h = H2mod(M,i)
s=size(M);
N=s(2);
L=lap(M);
range=1:N;
Li=L(range~=i,range~=i);
sigma = lyap(Li,-eye(N-1));
h = sqrt(trace(sigma));
