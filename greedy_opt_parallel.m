dellapool=parpool('local', str2num(getenv('PROCS'))) ;

numsigs_permove=1;
nummoves=1000;

N=20;
b=1;

timesteps=40;
evolution=zeros(N,timesteps);
groupspeed=zeros(1,timesteps);
groupconsensus=zeros(1,timesteps);
corrlengths=zeros(1,timesteps);

for t=1:timesteps

if t==1
    strategy=2*randi([1,(N-2)/2],1,N);
%     strategy=2*ones(1,N);
%     strategy(1)=10; 
    [probeaten, probgettoeat, meanlambda, meanH2, meancorrlength]=signalingevents_wgroupprops(strategy,numsigs_permove,nummoves,radius,b,T);
    groupspeed(1)=meanlambda;
    groupconsensus(1)=meanH2;
    corrlengths(1)=meancorrlength;
else strategy=evolution(:,t-1);
    strategy=reshape(strategy,1,N);

    M=makenet(strategy);

    v1=1./power(speeds(M,'both','max'),1-alpha);
    v2=1./power(vars(M,'additive','max'),alpha);
    perf=probeaten;

    newstrategy=strategy;

    for i=1:N
        choosestrat=zeros(1,3);
        choosestrat(2)=perf(i);

        if strategy(i)==N-1
            choosestrat(3)=perf(i);
        else 
            upone=strategy;
            upone(i)=strategy(i)+1;

            M=makenet(upone);
            v1=1./power(speeds(M,'both','max'),1-alpha);
            v2=1./power(vars(M,'additive','max'),alpha);
            choosestrat(3)=v1(i)*v2(i);
        end

        if strategy(i)==1
            choosestrat(1)=perf(i);
        else
            downone=strategy;
            downone(i)=strategy(i)-1;

            M=makenet(downone);
            v1=1./power(speeds(M,'both','max'),1-alpha);
            v2=1./power(vars(M,'additive','max'),alpha);
            choosestrat(1)=v1(i)*v2(i);   
        end 

        if choosestrat(1)>choosestrat(2)
            newstrategy(i)=strategy(i)-1;
        end
        if choosestrat(3)>choosestrat(2)
            newstrategy(i)=strategy(i)+1;
        end

    end
    M=makenet(newstrategy);
    L=lap(M);
    [~,d]=eig(L);
    d=diag(d);
    d=d(2:length(d));
    lambda=min(real(d));
    groupspeed(t)=1/((exp(-lambda))^(1-alpha))*1/((H2norm(M,'additive')/sqrt(N))^(alpha));
    evolution(:,t)=newstrategy;
end
end


% plot(1:timesteps,evolution(1,:))
% hold all
% for i=2:N
%     plot(1:timesteps,evolution(i,:))
% end
% hold off
% axis([1 timesteps 1 N-1])
% 
% figure
% plot(evolution(:,timesteps))
% xlabel('Node')
% ylabel('Strategy')
% figure
% plot(evolution(1,:))
% xlabel('Time')
% ylabel('Strategy')
% figure
% plot(groupperf)
% xlabel('Time')
% ylabel('Group Performance')
