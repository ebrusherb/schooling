numworkers=str2num(getenv('PROCS'));
dellacluster=parcluster('local');
dellacluster.JobStorageLocation=strcat('/scratch/network/brush/tmp/',getenv('SLURM_JOB_ID'));
dellapool=parpool(dellacluster, numworkers) ;

numsigs_permove=str2num(getenv('numsigs'));
nummoves=str2num(getenv('nummoves'));

N=20;
b=1;

strats=1:1:(N-1);
L=length(strats);

Tvals=[.01 .1 1 10];
Nt=length(Tvals);

radvals=[0:.1:1.4];
% radvals=.1;
Nr=length(radvals);

radius=radvals(2);
T=Tvals(3);
strategy=[1;19*ones(N-1,1)];

t0=tic;
[probeaten, probgettoeat]=signalingevents(strategy,numsigs_permove,nummoves,radius,b,T);
t=toc(t0);

t0=tic;
[probeatenshort, probgettoeatshort]=signalingevents_parallel(strategy,numsigs_permove,nummoves,radius,b,T);
tshort=toc(t0);

t0=tic;
[probeatenlong, probgettoeatlong]=signalingevents_parallel(strategy,numsigs_permove*N,nummoves,radius,b,T);
tlong=toc(t0);

filename=strcat('/home/brush/schooling_consensus/ESStest','.mat');
save(filename,'t','tshort','tlong','probeaten','probeatenshort','probeatenlong');

delete(dellapool);

exit ;
