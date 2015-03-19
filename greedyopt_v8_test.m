numworkers=str2num(getenv('PROCS')); %#ok<ST2NM>
dellacluster=parcluster('local');
dellacluster.JobStorageLocation=strcat('/scratch/network/brush/tmp/',getenv('SLURM_JOB_ID'));
dellapool=parpool(dellacluster, numworkers) ;

aid = getenv('SLURM_ARRAY_TASK_ID');
aid = str2num(aid); %#ok<ST2NM>

numsigs_permove=1;
nummoves=1000;

N=20;
b=1;
radvals=[.1 .2 .7];
% Nr=length(radvals);
T=1;

timesteps=200;
its=10;

evolution_eaten=zeros(its,N,timesteps);
groupconsensus_eaten=zeros(its,timesteps);
groupconsensusforced_eaten=zeros(its,timesteps);
corrlengths_eaten=zeros(its,timesteps);
corrlengthsforced_eaten=zeros(its,timesteps);
disconnected_eaten=zeros(its,timesteps);
changeorder_eaten=zeros(its,N,timesteps);

evolution_gettoeat=zeros(its,N,timesteps);
groupconsensus_gettoeat=zeros(its,timesteps);
groupconsensusforced_gettoeat=zeros(its,timesteps);
corrlengths_gettoeat=zeros(its,timesteps);
corrlengthsforced_gettoeat=zeros(its,timesteps);
disconnected_gettoeat=zeros(its,timesteps);

maxt_eaten=zeros(its,1);
maxt_gettoeat=zeros(its,1);

% for ir=1:Nr
    radius=radvals(aid);
for num=1:its
    init=randi([1 N/2-1],1,N);
    t=1;
    while t<=timesteps

    if t==1
        strategy=init;
    %     strategy=2*ones(1,N);
    %     strategy(1)=10; 
        evolution_eaten(num,:,t)=strategy;
        t=t+1;
    else
        strategy=evolution_eaten(num,:,t-1);
        strategy=reshape(strategy,1,N);

        [probeaten, ~]=signalingevents_parallel(strategy,numsigs_permove,nummoves,radius,b,T);
        perf=1-probeaten;
        newstrategy=strategy;

%         i=randsample(1:N,1);
%         for i=1:N
        changeorder=randsample(1:N,N,'false');
        changeorder_eaten(i,:,t)=changeorder;     
        for i=changeorder
            [probeaten, ~]=signalingevents_parallel(newstrategy,numsigs_permove,nummoves,radius,b,T);
            perf=1-probeaten;
            
            choosestrat=zeros(1,3);
            choosestrat(2)=perf(i);

            if strategy(i)==N-1
                choosestrat(3)=perf(i);
            else 
                upone=newstrategy;
                upone(i)=newstrategy(i)+1;

                [probeaten, ~]=signalingevents_parallel(upone,numsigs_permove,nummoves,radius,b,T);
                choosestrat(3)=1-probeaten(i);
            end

            if strategy(i)==1
                choosestrat(1)=perf(i);
            else
                downone=newstrategy;
                downone(i)=newstrategy(i)-1;

                [probeaten, ~]=signalingevents_parallel(downone,numsigs_permove,nummoves,radius,b,T);
                choosestrat(1)=1-probeaten(i);   
            end 

            if choosestrat(1)>choosestrat(2)
                newstrategy(i)=newstrategy(i)-1;
            end
            if choosestrat(3)>choosestrat(2)
                newstrategy(i)=newstrategy(i)+1;
            end

        end

        evolution_eaten(num,:,t)=newstrategy;
        if sum(evolution_eaten(num,:,t)==evolution_eaten(num,:,t-1))==N
            evolution_eaten(num,:,(t+1):end)=repmat(col(newstrategy),1,timesteps-t);
            t_eaten=t;
            t=timesteps+1;
        else
            t_eaten=t;
            t=t+1; 
        end
    end
    end
    maxt_eaten(num)=t_eaten;

    t=1;
    while t<=timesteps

    if t==1
        strategy=init;
    %     strategy=2*ones(1,N);
    %     strategy(1)=10; 
        evolution_gettoeat(num,:,t)=strategy;
        t=t+1;
    else
        strategy=evolution_gettoeat(num,:,t-1);
        strategy=reshape(strategy,1,N);

        [~, probgettoeat]=signalingevents_parallel(strategy,numsigs_permove,nummoves,radius,b,T);
        perf=probgettoeat;
        newstrategy=strategy;

%         i=randsample(1:N,1);
%         for i=1:N
        changeorder=randsample(1:N,N,'false');
        changeorder_eaten(i,:,t)=changeorder;     
        for i=changeorder
            choosestrat=zeros(1,3);
            choosestrat(2)=perf(i);

            if strategy(i)==N-1
                choosestrat(3)=perf(i);
            else 
                upone=strategy;
                upone(i)=newstrategy(i)+1;

                [~, probgettoeat]=signalingevents_parallel(upone,numsigs_permove,nummoves,radius,b,T);
                choosestrat(3)=probgettoeat(i);
            end

            if strategy(i)==1
                choosestrat(1)=perf(i);
            else
                downone=strategy;
                downone(i)=newstrategy(i)-1;

                [~, probgettoeat]=signalingevents_parallel(downone,numsigs_permove,nummoves,radius,b,T);
                choosestrat(1)=1-probgettoeat(i);   
            end 

            if choosestrat(1)>choosestrat(2)
                newstrategy(i)=newstrategy(i)-1;
            end
            if choosestrat(3)>choosestrat(2)
                newstrategy(i)=newstrategy(i)+1;
            end

        end
        evolution_gettoeat(num,:,t)=newstrategy;
        if sum(evolution_gettoeat(num,:,t)==evolution_gettoeat(num,:,t-1))==N
            evolution_gettoeat(num,:,(t+1):end)=repmat(col(newstrategy),1,timesteps-t);
            t_gettoeat=t;
            t=timesteps+1;
        else 
            t_gettoeat=t;
            t=t+1;
        end
    end
    end
    maxt_gettoeat(num)=t_gettoeat;
end

end

filename=strcat('/home/brush/schooling_consensus/greedyopt','_T=',num2str(T),'_nummoves=',num2str(nummoves),'_numpermove=',num2str(numsigs_permove),'_timesteps=',num2str(timesteps),'.mat');
save(filename,'evolution_eaten','evolution_gettoeat','disconnected_gettoeat','radvals','maxt_eaten','maxt_gettoeat')
'Greedy opt is done. Group props are not.' %#ok<NOPTS>

% for ir=1:Nr
for num=1:its
    for t=1:maxt_eaten(num)
        strategy_eaten=evolution_eaten(num,:,t);
        [meanH2, meanH2_forced, corrlength, corrlength_forced, disconnectedcount]=groupprops(strategy_eaten,numsigs_permove,nummoves,radius,b,T);
        groupconsensus_eaten(num,t)=meanH2;
        groupconsensusforced_eaten(num,t)=meanH2_forced;
        corrlengths_eaten(num,t)=corrlength;
        corrlengthsforced_eaten(num,t)=corrlength_forced;
        disconnected_eaten(num,t)=disconnectedcount;
    end
    groupconsensus_eaten(num,(t+1):end)=meanH2;
    groupconsensusforced_eaten(num,(t+1):end)=meanH2_forced;
    corrlengths_eaten(num,(t+1):end)=corrlength;
    corrlengthsforced_eaten(num,(t+1):end)=corrlength_forced;
    disconnected_eaten(num,(t+1):end)=disconnectedcount;
    for t=1:maxt_gettoeat(num)
        strategy_gettoeat=evolution_gettoeat(num,:,t);
        [meanH2, meanH2_forced, corrlength, corrlength_forced, disconnectedcount]=groupprops(strategy_gettoeat,numsigs_permove,nummoves,radius,b,T);
        groupconsensus_gettoeat(num,t)=meanH2;
        groupconsensusforced_gettoeat(num,t)=meanH2_forced;
        corrlengths_gettoeat(num,t)=corrlength;
        corrlengthsforced_gettoeat(num,t)=corrlength_forced;
        disconnected_gettoeat(num,t)=disconnectedcount;
    end
    groupconsensus_gettoeat(num,(t+1):end)=meanH2;
    groupconsensusforced_gettoeat(num,(t+1):end)=meanH2_forced;
    corrlengths_gettoeat(num,(t+1):end)=corrlength;
    corrlengthsforced_gettoeat(num,(t+1):end)=corrlength_forced;
    disconnected_gettoeat(num,(t+1):end)=disconnectedcount;
end
% end

filename=strcat('/home/brush/schooling_consensus/greedyopt','_T=',num2str(T),'_nummoves=',num2str(nummoves),'_numpermove=',num2str(numsigs_permove),'_rad=',num2str(radius),'_timesteps=',num2str(timesteps),'.mat');
save(filename,'evolution_eaten','groupconsensus_eaten','groupconsensusforced_eaten','corrlengths_eaten','corrlengthsforced_eaten','evolution_gettoeat','groupconsensus_gettoeat','groupconsensusforced_gettoeat','corrlengths_gettoeat','corrlengthsforced_gettoeat','disconnected_eaten','disconnected_gettoeat','radvals','maxt_eaten','maxt_gettoeat')
'Greedy opt is done. Group props are too.' %#ok<NOPTS>

delete(dellapool);

exit ;
