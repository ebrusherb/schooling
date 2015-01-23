numworkers=str2num(getenv('PROCS')); %#ok<*ST2NM>
dellacluster=parcluster('local');
dellacluster.JobStorageLocation=strcat('/scratch/network/brush/tmp/',getenv('SLURM_JOB_ID'));
dellapool=parpool(dellacluster, numworkers) ;

numsigs_permove=str2num(getenv('numsigs'));
nummoves=str2num(getenv('nummoves'));

aid = getenv('SLURM_ARRAY_TASK_ID');
aid = str2num(aid);

N=50;
b=1;

strats=10:1:(N-1);
L=length(strats);

Tvals=[1]; %#ok<*NBRAK>
Nt=length(Tvals);

radvals=[0:.1:1.4];
% radvals=.1;
Nr=length(radvals); %#ok<NASGU>

ESSseaten=cell(Nt);
ESSsgettoeat=cell(Nt);

storefitnesseaten=cell(Nt);
storerhoeaten=cell(Nt);

storefitnessgettoeat=cell(Nt);
storerhogettoeat=cell(Nt);

t0=tic;

for q=1:Nt
%     for p=1:Nr
        T=Tvals(q);
        radius=radvals(aid);
%         radius=radvals(p);

        fitnesseaten=zeros(L,L);
        rhoeaten=zeros(L,L);

        fitnessgettoeat=zeros(L,L);
        rhogettoeat=zeros(L,L);

        featen=zeros(L,L,N-1);
        geaten=zeros(L,L,N-1);

        fgettoeat=zeros(L,L,N-1);
        ggettoeat=zeros(L,L,N-1);
        parfor ind=1:(L*L*(N-1))
%         parfor ind=1:5
            [u,v,k]=ind2sub([L,L,N-1],ind);
            
            resident=strats(u); %#ok<PFBNS>
            invader=strats(v);

            strategy=resident*ones(1,N);
            strategy(N+1-(1:k))=invader;
            [probeaten, probgettoeat]=signalingevents(strategy,numsigs_permove,nummoves,radius,b,T);

            perfeaten=1-probeaten;
            featen(ind)=mean(perfeaten(N+1-(1:k)));
            geaten(ind)=mean(perfeaten(N+1-((k+1):N)));

            perfgettoeat=probgettoeat;
            fgettoeat(ind)=mean(perfgettoeat(N+1-(1:k)));
            ggettoeat(ind)=mean(perfgettoeat(N+1-((k+1):N)));

        end


        for ind=1:L*L
            [u,v]=ind2sub([L,L],ind);
            where=zeros(N-1,1);
            for k=1:N-1
                where(k)=sub2ind([L,L,N-1],u,v,k);
            end
            featen_now=featen(where);
            geaten_now=geaten(where);

            fgettoeat_now=fgettoeat(where);
            ggettoeat_now=ggettoeat(where);

            fitnesseaten(ind)=featen_now(1)/geaten_now(1);
            fitnessgettoeat(ind)=fgettoeat_now(1)/ggettoeat_now(1);

            ratio=geaten_now./featen_now;
            tosum=ones(1,N-1);
            for k=1:(N-1)
                tosum(k)=prod(ratio(1:k));
            end
            rhoeaten(ind)=1/(1+sum(tosum));

            ratio=ggettoeat_now./fgettoeat_now;
            tosum=ones(1,N-1);
            for k=1:(N-1)
                tosum(k)=prod(ratio(1:k));
            end
            rhogettoeat(ind)=1/(1+sum(tosum));
        end  
        
        storefitnesseaten{q}=fitnesseaten;
        storerhoeaten{q}=rhoeaten;
        
        storefitnessgettoeat{q}=fitnessgettoeat;
        storerhogettoeat{q}=rhogettoeat;

        eqstratseaten=eq_strats(N,fitnesseaten,rhoeaten);
        ESSseaten{q}=eqstratseaten{2};

        eqstratsgettoeat=eq_strats(N,fitnessgettoeat,rhogettoeat);
        ESSsgettoeat{q}=eqstratsgettoeat{2};
%     end
end

t=toc(t0);
disp(t);

filename=strcat('/home/brush/schooling_consensus/ESS_bigN','_nummoves=',num2str(nummoves),'_numpermove=',num2str(numsigs_permove),'_rad=',num2str(radius),'.mat');
save(filename,'strats','radvals','Tvals','storefitnesseaten','storerhoeaten','ESSseaten','storefitnessgettoeat','storerhogettoeat','ESSsgettoeat');

delete(dellapool);

exit ;
