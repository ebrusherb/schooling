
% dellapool=parpool('local', str2num(getenv('PROCS'))) ;
t0=tic;
numsigs_permove=1;
nummoves=100;


N=20;
b=1;

strats=2:2:(N-1);
L=length(strats);

Tvals=1;
Nt=length(Tvals);

radvals=[.0 .1 .5 1.5];
% radvals=.1;
Nr=length(radvals);

ESSseaten=cell(Nt,Nr);
ESSsgettoeat=cell(Nt,Nr);

q=1;
p=1;
        T=Tvals(q);
        radius=radvals(p);

        fitnesseaten=zeros(L,L);
        rhoeaten=zeros(L,L);

        fitnessgettoeat=zeros(L,L);
        rhogettoeat=zeros(L,L);

        featen=zeros(L,L,N-1);
        geaten=zeros(L,L,N-1);

        fgettoeat=zeros(L,L,N-1);
        ggettoeat=zeros(L,L,N-1);
%         for ind=1:(L*L*(N-1))
        for ind=1:2
            [u,v,k]=ind2sub([L,L,N-1],ind);
            
            resident=strats(u);
            invader=strats(v);

            strategy=resident*ones(1,N);
            strategy(1:k)=invader;
            [probeaten, probgettoeat]=signalingevents(strategy,numsigs_permove,nummoves,radius,b,T);

            perfeaten=1-probeaten;
            featen(ind)=mean(perfeaten(1:k));
            geaten(ind)=mean(perfeaten((k+1):N));

            perfgettoeat=probgettoeat;
            fgettoeat(ind)=mean(perfgettoeat(1:k));
            ggettoeat(ind)=mean(perfgettoeat((k+1):N));

        end

t=toc(t0);
% save('/home/brush/schooling_consensus/ESSoutput.mat','strats','radvals','fitnesseaten','rhoeaten','ESSseaten','fitnessgettoeat','rhogettoeat','ESSsgettoeat');

% delete(dellapool);

% exit ;
