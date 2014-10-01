function strat=dimorph(fitness)
s=size(fitness);
N=s(2);
strat=[];
count=0;
for i=1:(N-1)
    for j=(i+1):N
        if fitness(i,j)>fitness(i,i) && fitness(j,i)>fitness(j,j)
            count=count+1;
            strat(count,:)=[i j];
        end
    end            
end
