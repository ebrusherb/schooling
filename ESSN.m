function Strat=ESSN(fitness,rho,N)
s=size(fitness);
S=s(2);
d=diag(fitness);
range=1:S;
strat=[];
for k=1:S
    m=max(fitness(k,range~=k));
    if fitness(k,k)>m
       strat=[strat k]; 
    end
    if m==fitness(k,k)
        v=range(fitness(k,:)==m);
        if sum((fitness(v,k)>d(v))==0)==0
            strat=[strat k];
        end
    end            
end
if nargin>1
    Strat=[];
   for k=strat
       f=max(rho(k,range~=k));
        if f<=1/N
            Strat=[Strat k];
        end
   end
else Strat=strat;
end