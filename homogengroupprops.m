numworkers=str2num(getenv('PROCS'));
dellacluster=parcluster('local');
dellacluster.JobStorageLocation=strcat('/scratch/network/brush/tmp/',getenv('SLURM_JOB_ID'));
dellapool=parpool(dellacluster, numworkers) ;

numsigs_permove=1;
nummoves=1000;

N=20;
b=1;
T=1;

homogenstrats=6:1:(N-1);
Lh=length(homogenstrats);

radvals=0:.1:1.4;
Nr=length(radvals);

groupspeed=zeros(Lh,Nr);
groupconsensus=zeros(Lh,Nr);
corrlengths=zeros(Lh,Nr);
disconnected=zeros(Lh,Nr);

for i=1:Nr
    radius=radvals(i);
    parfor j=1:Lh
        strat=homogenstrats(j);
        strategy=strat*ones(N,1);
        [meanlambda, meanH2, meancorrlength, disconnectedcount]=groupprops(strategy,numsigs_permove,nummoves,radius,b,T);
        groupspeed(j,i)=meanlambda;
        groupconsensus(j,i)=meanH2;
        corrlengths(j,i)=meancorrlength;
        disconnected(j,i)=disconnectedcount;
    end
end

save('/home/brush/schooling_consensus/homogengroupprops.mat','groupspeed','groupconsensus','corrlengths','disconnected','homogenstrats','Lh','radvals','Nr','N')

delete(dellapool);

exit ;
