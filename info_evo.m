%% initialize positions and distance matrix
N=100;

positions=zeros(N,2);
positions(:,1)=rand(1,N);
positions(:,2)=rand(1,N);
d=squareform(pdist(positions));

%% set up parameters and vectors for evolution
generations=10;

dist_its=zeros([generations N N ]);
strategy=zeros([generations+1 N]);
%strategy(1,:)=randi(N-1,1,N);
strategy(1,1:12)=5;
strategy(1,13:N)=15;
H2vals=zeros(1,generations);

signal=0;
runs=100; 

%% evolve

for k=1:generations
    
    M=zeros(N);
    for i=1:N
        [ranked order]=sort(d(i,:));
        neighbors=order(2:strategy(k,i)+1);
        M(i,neighbors)=1;
    end
    L=lap(M);
    H2vals(1,k)=H2norm(L);   
    
    for chosen=1:N
        x=zeros([N runs]);
        x(:,1)=rand(1,N);
        x(chosen,1)=signal;
        for i=2:runs
            x(:,i)=dynamics(x(:,i-1),M);
            x(chosen,i)=signal;
            x(:,i)=x(:,i)+.1*randn(N,1);
        end    
        for i=1:N
            dist_its(k,chosen,i)=sqrt(sum(power(x(i,:)-x(chosen,:),2)));
        end
    end
    avg_dists=zeros(1,N);
    for i=1:N
        avg_dists(1,i)=mean2(dist_its(k,(1:N)~=i,i));
    end
    strategy(k+1,:)=nextgen(strategy(k,:),avg_dists);
end

%% trying to find an ESS strategy
    N=100;

    runs=100; 
    signal=0;
    
    positions=zeros(N,2);
    positions(:,1)=rand(1,N);
    positions(:,2)=rand(1,N);
    d=squareform(pdist(positions));
    
    
    
 resident=1:10:N;
 invader=1:10:N;
 s=size(resident);
 s=s(2);
 
 ESS=zeros(s,s);
 
 for count1=1:s
     for count2=1:s
    
    play1=resident(count1);
    play2=invader(count2);
%     strategy=randi(N-1,1,N);
%     strategy=99*ones(1,N);
%     strategy(1,1)=5;
    strategy=play1*ones(1,N);
    strategy(1,1)=play2;
    devs=zeros([N N runs]);

    M=zeros(N);
    for i=1:N
        [ranked order]=sort(d(i,:));
        neighbors=order(2:strategy(1,i)+1);
        M(i,neighbors)=1;
    end
    
    x=zeros([N N runs]);
    for chosen=1:N
        x(chosen,:,1)=rand(1,N)*2-1;
        x(chosen,chosen,1)=signal;
        for time=2:runs
            x(chosen,:,time)=dynamics(transpose(x(chosen,:,time-1)),M);
            x(chosen,chosen,time)=signal;
            x(chosen,:,time)=x(chosen,:,time)+.1*randn(1,N);
            devs(chosen,:,time)=x(chosen,:,time)-x(chosen,chosen,time);
        end    
    end
    avg_dists=zeros(1,N);
    var_dists=zeros(1,N);
    for i=1:N
        v=zeros(1,(N-1)*runs);
        count=0;
        for j=setdiff(1:N,i)
            v(1,(1:runs)+count*runs)=devs(j,i,:);
            count=count+1;
        end
        avg_dists(1,i)=mean2(reshape(devs(setdiff(1:N,i),i,:),1,[]));
        var_dists(1,i)=var(reshape(devs(setdiff(1:N,i),i,:),1,[]));
    end
    
%     scatter(strategy,avg_dists)
%     scatter(strategy,var_dists)
%     scatter(avg_dists,var_dists)
    ESS(count1,count2)=mean(avg_dists(2:N))-avg_dists(1);
     end
 end 
    
    
%% 
N=5;

q=zeros(N-1,N);
q(1,1:2)=[1/sqrt(2) -1/sqrt(2)];
for j=2:(N-1)
   v=rand(1,N);
   v=v-sum(v)/N;
   for k=1:(j-1)
      v=v-(v*transpose(q(k,:)))*q(k,:); 
   end
   v=v/norm(v,2);
   q(j,:)=v;
end

    its=500;
    
    H2vals=zeros(1,its);
    H2valsmod=zeros(1,its);
    min_rates=zeros(1,its);
    max_rates=zeros(N,its);
    
    for k=1:its
    
    S=2;
    
    while S>1
        
    positions=zeros(N,2);
    positions(:,1)=rand(1,N);
    positions(:,2)=rand(1,N);
    d=squareform(pdist(positions));
    
    strategy=randi(N-1,1,N);

    M=zeros(N);
    for i=1:N
        [ranked order]=sort(d(i,:));
        neighbors=order(2:strategy(1,i)+1);
        M(i,neighbors)=1;
    end
    
    [S,~]=graphconncomp(sparse(M));
    
    end
    
    H2vals(k)=H2norm(M);
    for i=1:N
    H2valsmod(i,k)=H2mod(M,i);
    end
    min_rates(k)=min(convergence2(d,strategy));
    max_rates(k)=max(convergence2(d,strategy));
    
    
    end
    
    scatter(H2vals,min_rates);
    scatter(mean(H2valsmod),H2vals)
    
 %% compare speed and variance trade off

 N=50;
 
 S=2;
    
    while S>1
        
    positions=zeros(N,2);
    positions(:,1)=rand(1,N);
    positions(:,2)=rand(1,N);
    d=squareform(pdist(positions));
    
    strategy=randi(N-1,1,N);

    M=zeros(N);
    for i=1:N
        [ranked order]=sort(d(i,:));
        neighbors=order(2:strategy(1,i)+1);
        M(i,neighbors)=1;
    end
    
    [S,~]=graphconncomp(sparse(M));
    
    end  
    
    Vars=vars(M);
    Speeds=speeds(M);
    scatter(Vars,Speeds)