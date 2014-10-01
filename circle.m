function M=circle(n,k)
M=zeros(n,n);
l=ceil(k/2);
indices=[1:n 1:n 1:n];
if k==1
    for i=1:(n-1)
        M(i,i+1)=1;
    end
    M(n,1)=1;
else 
    for i=1:n
        M(i,indices(i+(1:l)+n))=1;
        M(i,indices(i-(1:(k-l))+n))=1;
    end
end

       