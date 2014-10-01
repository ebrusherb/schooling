function strats=eq_strats(N,fitness,rho)
%input fitness / rho have residents in rows and invaders in columns
s=size(fitness,2);
d=diag(fitness);
range=1:s;
ESS=[];
NIS=[];
branch=[];
if nargin<3
    rho=1/N*ones(N,N);
end

superf=fitness(1:(s-1),2:s);
superf=diag(superf);
superf=[1; superf];
superr=rho(1:(s-1),2:s);
superr=diag(superr);
superr=[1/N; superr];
increase=(superf>=1-0.0001).*(superr>=1/N-0.0001); %increase(i)=1 if you'll increase TO i

subf=fitness(2:s,1:(s-1));
subf=diag(subf);
subf=[subf; 1];
subr=rho(2:s,1:(s-1));
subr=diag(subr);
subr=[subr; 1/N];
decrease=(subf>=1-0.0001).*(subr>=1/N-0.00001); %decrease(i)=1 if you'll decrease TO i

converge=range(increase.*decrease==1); %converge(i)=1 if you'll both increase and decrease TO i

% clear superf supperr subf subr 

for m=converge
    if m>1 && m<s
        if increase(m+1)==1 && decrease(m-1)==1
            branch=[branch m];
        end
    end
end

for k=1:s
    
    m=max(fitness(k,range~=k)); %the biggest fitness of invaders to k
    f=max(rho(k,range~=k)); %the biggest fixation probability of invaders to k
    if fitness(k,k)>m && f<=1/N+0.0001 %if k's fitness is bigger than its invaders and no one can fixate 
       ESS=[ESS k]; %k is an ESS
    end
    if m==fitness(k,k) %if there is an invader that can do as well
        v=range(fitness(k,:)==m);
        if sum((fitness(v,k)>d(v))==0)==0 && min(rho(v,k))>=1/N-0.0001 %then check if k as an invader can do better against those invaders
            ESS=[ESS k];
        end
    end
    
    f=min(rho(range~=k,k)); %the lowest fixation probability of k as an invader
    if sum((d(range~=k)<fitness(range~=k,k))==1)==(s-1) && f>=1/N-0.0001 %if k as an invader does better than all residents and it can come to fixation
       NIS=[NIS k]; %k is an NIS
    end
    if sum((d(range~=k)==fitness(range~=k,k))==1)>=1 %if there is a resident k does equally well against
        v=range((d(range~=k)==fitness(range~=k,k))); 
        f=min(rho(v,k)); %minimum fixation probability against those residents
        if f>=1/N-0.0001 && max(fitness(k,v))<d(k)+0.0001 && max(rho(k,v))<=1/N+0.0001 %then check if k as a resident can do better against those residents
            NIS=[NIS k];
        end
    end
end

% for k=2:(s-1)
%     
%     if (superf(k+1)>=1-0.0001).*(superr(k+1)>=1/N-0.0001)==0 && (subf(k-1)>=1-0.0001).*(subr(k-1)>=1/N-0.00001)==0
%         ESS=[ESS k];
%     end

CSS=intersect(ESS,NIS);

strats=cell(1,5);
strats{1}=converge;
strats{2}=ESS;
strats{3}=NIS;
strats{4}=CSS;
strats{5}=branch;

