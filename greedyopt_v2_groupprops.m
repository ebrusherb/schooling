numworkers=str2num(getenv('PROCS'));
dellacluster=parcluster('local');
dellacluster.JobStorageLocation=strcat('/scratch/network/brush/tmp/',getenv('SLURM_JOB_ID'));
dellapool=parpool(dellacluster, numworkers) ;

load /home/brush/schooling_consensus/greedyopt.mat

numsigs_permove=1;
nummoves=1000;

N=20;
b=1;
radius=.5;
T=1;

timesteps=50;
its=50;

groupspeed_eaten=zeros(its,timesteps);
groupconsensus_eaten=zeros(its,timesteps);
corrlengths_eaten=zeros(its,timesteps);
disconnected_eaten=zeros(its,timesteps);

groupspeed_gettoeat=zeros(its,timesteps);
groupconsensus_gettoeat=zeros(its,timesteps);
corrlengths_gettoeat=zeros(its,timesteps);
disconnected_gettoeat=zeros(its,timesteps);

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
end

save('/home/brush/schooling_consensus/greedyopt_groupprops.mat','groupspeed_eaten','groupconsensus_eaten','corrlengths_eaten','groupspeed_gettoeat','groupconsensus_gettoeat','corrlengths_gettoeat','disconnected_eaten','disconnected_gettoeat')

delete(dellapool);

exit ;
