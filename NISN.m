function Strat=NISN(fitness,rho,N)
s=size(fitness);
S=s(2);
d=diag(fitness);
range=1:S;
strat=[];
for k=1:S
    if sum((d(range~=k)<fitness(range~=k,k))==1)==(S-1)
       strat=[strat k]; 
    end
    if sum((d(range~=k)==fitness(range~=k,k))==1)>=1
        v=range((d(range~=k)==fitness(range~=k,k)));
        if max(fitness(k,v))<d(k)
            strat=[strat k];
        end
    end            
end
if nargin>1
    Strat=[];
   for k=strat
       f=max(rho(range~=k,k));
        if f>=1/N
            Strat=[Strat k];
        end
   end
else Strat=strat;
end