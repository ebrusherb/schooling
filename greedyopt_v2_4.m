numworkers=str2num(getenv('PROCS')); %#ok<ST2NM>
dellacluster=parcluster('local');
dellacluster.JobStorageLocation=strcat('/scratch/network/brush/tmp/',getenv('SLURM_JOB_ID'));
dellapool=parpool(dellacluster, numworkers) ;

numsigs_permove=1000;
nummoves=1;

N=20;
b=1;
radius=.1;
T=1;

timesteps=50;
its=20;

evolution_eaten=zeros(its,N,timesteps);
groupspeed_eaten=zeros(its,timesteps);
groupconsensus_eaten=zeros(its,timesteps);
corrlengths_eaten=zeros(its,timesteps);
disconnected_eaten=zeros(its,timesteps);

evolution_gettoeat=zeros(its,N,timesteps);
groupspeed_gettoeat=zeros(its,timesteps);
groupconsensus_gettoeat=zeros(its,timesteps);
corrlengths_gettoeat=zeros(its,timesteps);
disconnected_gettoeat=zeros(its,timesteps);

for num=1:its

for t=1:timesteps
    init=randi([1 N-1],1,N);
if t==1
    strategy=init;
%     strategy=2*ones(1,N);
%     strategy(1)=10; 
    evolution_eaten(num,:,t)=strategy;
else
    strategy=evolution_eaten(num,:,t-1);
    strategy=reshape(strategy,1,N);

    [probeaten, ~]=signalingevents_parallel(strategy,numsigs_permove,nummoves,radius,b,T);
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

            [probeaten, ~]=signalingevents_parallel(upone,numsigs_permove,nummoves,radius,b,T);
            choosestrat(3)=1-probeaten(i);
        end

        if strategy(i)==1
            choosestrat(1)=perf(i);
        else
            downone=strategy;
            downone(i)=strategy(i)-1;

            [probeaten, ~]=signalingevents_parallel(downone,numsigs_permove,nummoves,radius,b,T);
            choosestrat(1)=1-probeaten(i);   
        end 

        if choosestrat(1)>choosestrat(2)
            newstrategy(i)=strategy(i)-1;
        end
        if choosestrat(3)>choosestrat(2)
            newstrategy(i)=strategy(i)+1;
        end

    end
    evolution_eaten(num,:,t)=newstrategy;
end

if t==1
    strategy=init;
%     strategy=2*ones(1,N);
%     strategy(1)=10; 
    evolution_gettoeat(num,:,t)=strategy;
else
    strategy=evolution_gettoeat(num,:,t-1);
    strategy=reshape(strategy,1,N);

    [~, probgettoeat]=signalingevents_parallel(strategy,numsigs_permove,nummoves,radius,b,T);
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

            [~, probgettoeat]=signalingevents_parallel(upone,numsigs_permove,nummoves,radius,b,T);
            choosestrat(3)=probgettoeat(i);
        end

        if strategy(i)==1
            choosestrat(1)=perf(i);
        else
            downone=strategy;
            downone(i)=strategy(i)-1;

            [~, probgettoeat]=signalingevents_parallel(downone,numsigs_permove,nummoves,radius,b,T);
            choosestrat(1)=1-probgettoeat(i);   
        end 

        if choosestrat(1)>choosestrat(2)
            newstrategy(i)=strategy(i)-1;
        end
        if choosestrat(3)>choosestrat(2)
            newstrategy(i)=strategy(i)+1;
        end

    end
    evolution_gettoeat(num,:,t)=newstrategy;
end
end
end

for num=1:its
    for t=1:timesteps
        strategy_eaten=evolution_eaten(num,:,t);
        [meanlambda, meanH2, meancorrlength, disconnectedcount]=groupprops_parallel(strategy_eaten,numsigs_permove,nummoves,radius,b,T);
        groupspeed_eaten(num,t)=meanlambda;
        groupconsensus_eaten(num,t)=meanH2;
        corrlengths_eaten(num,t)=meancorrlength;
        disconnected_eaten(num,t)=disconnectedcount;
        strategy_gettoeat=evolution_gettoeat(num,:,t);
        [meanlambda, meanH2, meancorrlength, disconnectedcount]=groupprops_parallel(strategy_gettoeat,numsigs_permove,nummoves,radius,b,T);
        groupspeed_gettoeat(num,t)=meanlambda;
        groupconsensus_gettoeat(num,t)=meanH2;
        corrlengths_gettoeat(num,t)=meancorrlength;
        disconnected_gettoeat(num,t)=disconnectedcount;
    end
% end

filename=strcat('/home/brush/schooling_consensus/greedyopt','_rad=',num2str(radius),'_T=',num2str(T),'_nummoves=',num2str(nummoves),'_numpermove=',num2str(numsigs_permove),'.mat');
save(filename,'evolution_eaten','groupspeed_eaten','groupconsensus_eaten','corrlengths_eaten','evolution_gettoeat','groupspeed_gettoeat','groupconsensus_gettoeat','corrlengths_gettoeat','disconnected_eaten','disconnected_gettoeat')

delete(dellapool);

exit ;
