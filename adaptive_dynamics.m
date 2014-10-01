N=20;

alpha=.5;
timesteps=30;
evolution=zeros(N,timesteps);
groupperf=zeros(1,timesteps);

for t=1:timesteps

if t==1
    %strategy=randi([2,N-1],1,N);
    strategy=2*ones(1,N);
    strategy(1)=10;
    evolution(:,t)=strategy;
    M=makenet(strategy);
    L=lap(M);
    [~,d]=eig(L);
    lambda=max(diag(d));
    groupperf(1)=1/((exp(-lambda))^(1-alpha))*1/((H2norm(M,'additive')/sqrt(N))^(alpha));
else strategy=evolution(:,t-1);
    strategy=reshape(strategy,1,N);


M=makenet(strategy);

v1=1./power(speeds(M,'both','max'),1-alpha);
v2=1./power(vars(M,'additive','max'),alpha);
perf=v1.*v2;

die=randi([1 N],1,1);

probs=perf/sum(perf);
draw=rand(1);

newstrategy=strategy;

M=makenet(newstrategy);
    L=lap(M);
    [~,d]=eig(L);
    d=diag(d);
    lambda=d(2);
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

% t=1;
% Mt=makenet(evolution(:,t));
% plot(t,H2norm(Mt),'o')
% hold all
% for t=2:timesteps
%     Mt=makenet(evolution(:,t));
%     plot(t,H2norm(Mt),'o')
% end
% hold off
