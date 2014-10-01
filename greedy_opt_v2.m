% dellapool=parpool('local', str2num(getenv('PROCS'))) ;
% matlabpool open local 2

numsigs_permove=1;
nummoves=1000;

N=20;
b=1;
radius=.5;
T=1;

timesteps=200;
its=10;

evolution_eaten=zeros(its,N,timesteps);
groupspeed_eaten=zeros(its,timesteps);
groupconsensus_eaten=zeros(its,timesteps);
corrlengths_eaten=zeros(its,timesteps);

for num=1:its

for t=1:timesteps

if t==1
    strategy=randi([1 N-1],1,N);
%     strategy=2*ones(1,N);
%     strategy(1)=10; 
    [~, ~, meanlambda, meanH2, meancorrlength]=signalingevents_wgroupprops_parallel(strategy,numsigs_permove,nummoves,radius,b,T);
    groupspeed_eaten(num,1)=meanlambda;
    groupconsensus_eaten(num,1)=meanH2;
    corrlengths_eaten(num,1)=meancorrlength;
    evolution_eaten(num,:,t)=strategy;
else
    strategy=evolution_eaten(num,:,t-1);
    strategy=reshape(strategy,1,N);

    [probeaten, ~, ~, ~, ~]=signalingevents_wgroupprops_parallel(strategy,numsigs_permove,nummoves,radius,b,T);
    perf=1-probeaten;
    newstrategy=strategy;

    for i=1:N
        choosestrat=zeros(1,3);
        choosestrat(2)=perf(i);

        if strategy(i)==N-1
            choosestrat(3)=perf(i);
        else 
            upone=strategy;
            upone(i)=strategy(i)+1;

            [probeaten, ~, ~, ~, ~]=signalingevents_wgroupprops_parallel(upone,numsigs_permove,nummoves,radius,b,T);
            choosestrat(3)=1-probeaten(i);
        end

        if strategy(i)==1
            choosestrat(1)=perf(i);
        else
            downone=strategy;
            downone(i)=strategy(i)-1;

            [probeaten, ~, ~, ~, ~]=signalingevents_wgroupprops_parallel(downone,numsigs_permove,nummoves,radius,b,T);
            choosestrat(1)=1-probeaten(i);   
        end 

        if choosestrat(1)>choosestrat(2)
            newstrategy(i)=strategy(i)-1;
        end
        if choosestrat(3)>choosestrat(2)
            newstrategy(i)=strategy(i)+1;
        end

    end
    [~, ~, meanlambda, meanH2, meancorrlength]=signalingevents_wgroupprops_parallel(newstrategy,numsigs_permove,nummoves,radius,b,T);
    groupspeed_eaten(num,t)=meanlambda;
    groupconsensus_eaten(num,t)=meanH2;
    corrlengths_eaten(num,t)=meancorrlength;
    evolution_eaten(num,:,t)=newstrategy;
end
end
end

evolution_gettoeat=zeros(its,N,timesteps);
groupspeed_gettoeat=zeros(its,timesteps);
groupconsensus_gettoeat=zeros(its,timesteps);
corrlengths_gettoeat=zeros(its,timesteps);

for num=1:its

for t=1:timesteps

if t==1
    strategy=randi([1 N-1],1,N);
%     strategy=2*ones(1,N);
%     strategy(1)=10; 
    [~, ~, meanlambda, meanH2, meancorrlength]=signalingevents_wgroupprops_parallel(strategy,numsigs_permove,nummoves,radius,b,T);
    groupspeed_gettoeat(num,1)=meanlambda;
    groupconsensus_gettoeat(num,1)=meanH2;
    corrlengths_gettoeat(num,1)=meancorrlength;
    evolution_gettoeat(num,:,t)=strategy;
else
    strategy=evolution_gettoeat(num,:,t-1);
    strategy=reshape(strategy,1,N);

    [~, probgettoeat, ~, ~, ~]=signalingevents_wgroupprops_parallel(strategy,numsigs_permove,nummoves,radius,b,T);
    perf=probgettoeat;
    newstrategy=strategy;

    for i=1:N
        choosestrat=zeros(1,3);
        choosestrat(2)=perf(i);

        if strategy(i)==N-1
            choosestrat(3)=perf(i);
        else 
            upone=strategy;
            upone(i)=strategy(i)+1;

            [~, probgettoeat, ~, ~, ~]=signalingevents_wgroupprops_parallel(upone,numsigs_permove,nummoves,radius,b,T);
            choosestrat(3)=probgettoeat(i);
        end

        if strategy(i)==1
            choosestrat(1)=perf(i);
        else
            downone=strategy;
            downone(i)=strategy(i)-1;

            [~, probgettoeat, ~, ~, ~]=signalingevents_wgroupprops_parallel(downone,numsigs_permove,nummoves,radius,b,T);
            choosestrat(1)=1-probgettoeat(i);   
        end 

        if choosestrat(1)>choosestrat(2)
            newstrategy(i)=strategy(i)-1;
        end
        if choosestrat(3)>choosestrat(2)
            newstrategy(i)=strategy(i)+1;
        end

    end
    [~, ~, meanlambda, meanH2, meancorrlength]=signalingevents_wgroupprops_parallel(newstrategy,numsigs_permove,nummoves,radius,b,T);
    groupspeed_gettoeat(num,t)=meanlambda;
    groupconsensus_gettoeat(num,t)=meanH2;
    corrlengths_gettoeat(num,t)=meancorrlength;
    evolution_gettoeat(num,:,t)=newstrategy;
end
end
end
% matlabpool close
save('/home/brush/schooling_consensus/greedyopt.mat','evolution_eaten','groupspeed_eaten','groupconsensus_eaten','corrlengths_eaten','evolution_gettoeat','groupspeed_gettoeat','groupconsensus_gettoeat','corrlengths_gettoeat')

delete(dellapool);

exit ;
