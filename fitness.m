function m=fitness(N,alpha)
strats=2:2:(N-1);
    L=length(strats);
    fitness=zeros(L,L);
    rho=zeros(L,L);
    

        for i=1:L
            for j=1:L
            resident=strats(i);
            invader=strats(j);

            f=zeros(1,N-1);
            g=zeros(1,N-1);

            for k=1:(N-1)
                strategy=resident*ones(1,N);
                strategy(1:k)=invader;
                M=makenet(strategy);
                v1=speeds(M,'both','mean');
                v2=vars(M,'additive','mean');
                perf=1./((1-alpha)*v1+alpha*v2);
                f(k)=mean(perf(1:k));
                g(k)=mean(perf((k+1):N));
                if k==1
                    fitness(i,j)=f(k)/g(k);
                end
        
            end
            
            ratio=g./f;
            tosum=ones(1,N-1);
            for k=1:(N-1)
                tosum(k)=prod(ratio(1:k));
            end
            rho(i,j)=1/(1+sum(tosum));
            
            end
        
        end
m{1}=fitness;
m{2}=rho;