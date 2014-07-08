N=20;

alpha=.5;
timesteps=40;
evolution=zeros(N,timesteps);
groupperf=zeros(1,timesteps);

for t=1:timesteps

if t==1
    %strategy=2*randi([1,(N-2)/2],1,N);
    strategy=2*ones(1,N);
    strategy(1)=10;
    evolution(:,t)=strategy;
    M=makenet(strategy);
    L=lap(M);
    [~,d]=eig(L);
    d=diag(d);
    d=d(2:length(d));
    lambda=min(real(d));
    groupperf(1)=1/((exp(-lambda))^(1-alpha))*1/((H2norm(M,'additive')/sqrt(N))^(alpha));
else strategy=evolution(:,t-1);
    strategy=reshape(strategy,1,N);


M=makenet(strategy);

v1=1./power(speeds(M,'both','max'),1-alpha);
v2=1./power(vars(M,'additive','max'),alpha);
perf=v1.*v2;

newstrategy=strategy;

for i=1:N
    choosestrat=zeros(1,3);
    choosestrat(2)=perf(i);
    
    if strategy(i)==N-2
        choosestrat(3)=perf(i);
    else 
    upone=strategy;
    upone(i)=strategy(i)+2;
    
    M=makenet(upone);
    v1=1./power(speeds(M,'both','max'),1-alpha);
    v2=1./power(vars(M,'additive','max'),alpha);
    choosestrat(3)=v1(i)*v2(i);
    end
    
    if strategy(i)==2
        choosestrat(1)=perf(i);
    else

    downone=strategy;
    downone(i)=strategy(i)-2;
 
    M=makenet(downone);
    v1=1./power(speeds(M,'both','max'),1-alpha);
    v2=1./power(vars(M,'additive','max'),alpha);
    choosestrat(1)=v1(i)*v2(i);   
    end 
   
    if choosestrat(1)>choosestrat(2)
        newstrategy(i)=strategy(i)-2;
    end
    if choosestrat(3)>choosestrat(2)
        newstrategy(i)=strategy(i)+2;
    end
    
end
M=makenet(newstrategy);
    L=lap(M);
    [~,d]=eig(L);
    d=diag(d);
    d=d(2:length(d));
    lambda=min(real(d));
    groupperf(t)=1/((exp(-lambda))^(1-alpha))*1/((H2norm(M,'additive')/sqrt(N))^(alpha));
evolution(:,t)=newstrategy;
end
end


plot(1:timesteps,evolution(1,:))
hold all
for i=2:N
    plot(1:timesteps,evolution(i,:))
end
hold off
axis([1 timesteps 1 N-1])

figure
plot(evolution(:,timesteps))
xlabel('Node')
ylabel('Strategy')
figure
plot(evolution(1,:))
xlabel('Time')
ylabel('Strategy')
figure
plot(groupperf)
xlabel('Time')
ylabel('Group Performance')
