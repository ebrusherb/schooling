function t=parallel_test(n,parallel)
t=zeros(1,n);
if strcmp(parallel,'yes')
    parfor i=1:n
        t(i)=i^2;
    end
else 
    for i=1:n
        t(i)=i^2;
    end
end