function v=vars(M,additiveoruniform,meanormax)
s=size(M);
N=s(2);
L=lap(M);
if nargin==1
    additiveoruniform='additive';
    meanormax='mean';
end
if strcmp(additiveoruniform,'additive')==1
noise=diag(sum(transpose(M)));
else
noise=eye(N);
end
V=zeros(N);
range=1:N;
if var(sum(transpose(M)))==0 %if everyone's the same don't bother looping over all receivers
    Li=L(range~=1,range~=1);
    sigma = lyap(Li,-noise(range~=1,range~=1));
    V = sqrt(diag(sigma));
    if strcmp(meanormax,'mean')==1
        v=mean(V)*ones(N,1);
    else
        v=max(V)*ones(N,1);
    end
end
if var(sum(transpose(M)))~=0
    for i=1:N
        Li=L(range~=i,range~=i);
        sigma = lyap(Li,-noise(range~=i,range~=i));
        V(range~=i,i) = sqrt(diag(sigma));
    end
    if strcmp(meanormax,'mean')==1
        v=mean(V,2);
    else 
        v=max(transpose(V));
        v=transpose(v);
    end
end