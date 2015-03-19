% numworkers=str2num(getenv('PROCS')); %#ok<ST2NM>
% dellacluster=parcluster('local');
% dellacluster.JobStorageLocation=strcat('/scratch/network/brush/tmp/',getenv('SLURM_JOB_ID'));
% dellapool=parpool(dellacluster, numworkers) ;

numsigs_permove=1;
nummoves=100;

N=20;
b=1;
T=1;

% homogenstrats=6:1:(N-1);
homogenstrats=6;
Lh=length(homogenstrats);

radvals=0:.1:1.4;
% radvals=[0 1.4];
Nr=length(radvals);

groupconsensus=zeros(Lh,Nr);
groupconsensus_forced=zeros(Lh,Nr);
corrlengths=zeros(Lh,Nr);
corrlengths_forced=zeros(Lh,Nr);
disconnected=zeros(Lh,Nr);

for i=1:Nr
    radius=radvals(i);
%     parfor j=1:Lh
    for j=1:Lh
        strat=homogenstrats(j);
        strategy=strat*ones(N,1);
        [meanH2, meanH2_forced, corrlength, corrlength_forced, disconnectedcount]=groupprops(strategy,numsigs_permove,nummoves,radius,b,T);
        groupconsensus(j,i)=meanH2;
        groupconsensus_forced(j,i)=meanH2_forced;
        corrlengths(j,i)=corrlength;
        corrlengths_forced(j,i)=corrlength_forced;
        disconnected(j,i)=disconnectedcount;
    end
end

% save('/home/brush/schooling_consensus/homogengroupprops.mat','groupconsensus','groupconsensus_forced','corrlengths','corrlengths_forced','disconnected','homogenstrats','Lh','radvals','Nr','N')
% 
% delete(dellapool);
% 
% exit ;
