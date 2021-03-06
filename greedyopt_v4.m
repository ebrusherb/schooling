numworkers=str2num(getenv('PROCS')); %#ok<ST2NM>
dellacluster=parcluster('local');
dellacluster.JobStorageLocation=strcat('/scratch/network/brush/tmp/',getenv('SLURM_JOB_ID'));
dellapool=parpool(dellacluster, numworkers) ;

numsigs_permove=1000;
nummoves=1;

N=20;
b=1;
radvals=[.1 .7];
Nr=length(radvals);
T=1;

timesteps=100;
its=20;

evolution_eaten=zeros(Nr,its,N,timesteps);
groupspeed_eaten=zeros(Nr,its,timesteps);
groupconsensus_eaten=zeros(Nr,its,timesteps);
corrlengths_eaten=zeros(Nr,its,timesteps);
disconnected_eaten=zeros(Nr,its,timesteps);

evolution_gettoeat=zeros(Nr,its,N,timesteps);
groupspeed_gettoeat=zeros(Nr,its,timesteps);
groupconsensus_gettoeat=zeros(Nr,its,timesteps);
corrlengths_gettoeat=zeros(its,timesteps);
disconnected_gettoeat=zeros(Nr,its,timesteps);

for ir=1:Nr
    radius=radvals(ir);
for num=1:its
t=1;
while t<=timesteps
    init=randi([1 N-1],1,N);
if t==1
    strategy=init;
%     strategy=2*ones(1,N);
%     strategy(1)=10; 
    evolution_eaten(ir,num,:,t)=strategy;
else
    strategy=evolution_eaten(ir,num,:,t-1);
    strategy=reshape(ir,strategy,1,N);

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
    evolution_eaten(ir,num,:,t)=newstrategy;
    if sum(evolution_eaten(ir,num,:,t)==evolution_eaten(ir,num,:,t-1)==N)
        evolution_eaten(ir,num,:,(t+1):end)=repmat(col(newstrategy),1,timesteps-t);
        t=timesteps+1;
    else 
        t=t+1;
    end
end
end

while t<=timesteps
t=1;
if t==1
    strategy=init;
%     strategy=2*ones(1,N);
%     strategy(1)=10; 
    evolution_gettoeat(ir,num,:,t)=strategy;
else
    strategy=evolution_gettoeat(ir,num,:,t-1);
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
    evolution_gettoeat(ir,num,:,t)=newstrategy;
    if sum(evolution_gettoeat(ir,num,:,t)==evolution_gettoeat(ir,num,:,t-1))==N
        evolution_gettoeat(ir,num,:,(t+1):end)=repmat(col(newstrategy),1,timesteps-t);
        t=timesteps+1;
    else 
        t=t+1;
    end
end
end
end
end

filename=strcat('/home/brush/schooling_consensus/greedyopt','_T=',num2str(T),'_nummoves=',num2str(nummoves),'_numpermove=',num2str(numsigs_permove),'.mat');
save(filename,'evolution_eaten','groupspeed_eaten','groupconsensus_eaten','corrlengths_eaten','evolution_gettoeat','groupspeed_gettoeat','groupconsensus_gettoeat','corrlengths_gettoeat','disconnected_eaten','disconnected_gettoeat','radvals')
'Greedy opt is done. Group props are not.'

for ir=1:Nr
for num=1:its
    for t=1:timesteps
        strategy_eaten=evolution_eaten(ir,num,:,t);
        [meanlambda, meanH2, meancorrlength, disconnectedcount]=groupprops(strategy_eaten,numsigs_permove,nummoves,radius,b,T);
        groupspeed_eaten(ir,num,t)=meanlambda;
        groupconsensus_eaten(ir,num,t)=meanH2;
        corrlengths_eaten(ir,num,t)=meancorrlength;
        disconnected_eaten(ir,num,t)=disconnectedcount;
        strategy_gettoeat=evolution_gettoeat(ir,num,:,t);
        [meanlambda, meanH2, meancorrlength, disconnectedcount]=groupprops(strategy_gettoeat,numsigs_permove,nummoves,radius,b,T);
        groupspeed_gettoeat(ir,num,t)=meanlambda;
        groupconsensus_gettoeat(ir,num,t)=meanH2;
        corrlengths_gettoeat(ir,num,t)=meancorrlength;
        disconnected_gettoeat(ir,num,t)=disconnectedcount;
    end
end
end

filename=strcat('/home/brush/schooling_consensus/greedyopt','_T=',num2str(T),'_nummoves=',num2str(nummoves),'_numpermove=',num2str(numsigs_permove),'.mat');
save(filename,'evolution_eaten','groupspeed_eaten','groupconsensus_eaten','corrlengths_eaten','evolution_gettoeat','groupspeed_gettoeat','groupconsensus_gettoeat','corrlengths_gettoeat','disconnected_eaten','disconnected_gettoeat','radvals')
'Greedy opt is done. Group props are too.'

delete(dellapool);

exit ;
