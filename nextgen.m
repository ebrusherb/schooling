function next=nextgen(strategy,fitness)
n=size(fitness,2);
l=floor(n/2);
[~,fitorder]=sort(fitness);
replace=fitorder(n-l+1:n);
next=strategy;
next(replace)=strategy(fitorder(1:l));

