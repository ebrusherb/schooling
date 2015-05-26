numworkers=str2num(getenv('PROCS')); %#ok<ST2NM>
dellacluster=parcluster('local');
dellacluster.JobStorageLocation=strcat('/scratch/network/brush/tmp/',getenv('SLURM_JOB_ID'));
dellapool=parpool(dellacluster, numworkers) ;

numsigs_permove=str2num(getenv('numsigs')); %#ok<ST2NM>
nummoves=str2num(getenv('nummoves'));%#ok<ST2NM>

N=40;
b=1;

strats=1:1:(N-1);
L=length(strats);

Tvals=1;
Nt=length(Tvals);

mvals=[1:2:N]; %#ok<NBRAK>
% radvals=.1;
Nm=length(mvals);

ESSseaten=cell(Nt,Nm);
ESSsgettoeat=cell(Nt,Nm);
ESSsboth=cell(Nt,Nm);
% ESSsgenerous=cell(Nt,Nm);

storefitnesseaten=cell(Nt,Nm);
storerhoeaten=cell(Nt,Nm);

storefitnessgettoeat=cell(Nt,Nm);
storerhogettoeat=cell(Nt,Nm);

storefitnessboth=cell(Nt,Nm);
storerhoboth=cell(Nt,Nm);

% storefitnessgenerous=cell(Nt,Nm);
% storerhogenerous=cell(Nt,Nm);

t0=tic;

for q=1:Nt
    for p=1:Nm
%         T=Tvals(q);
        m=mvals(p);

        fitnesseaten=zeros(L,L);
        rhoeaten=zeros(L,L);

        fitnessgettoeat=zeros(L,L);
        rhogettoeat=zeros(L,L);
        
        fitnessboth=zeros(L,L);
        rhoboth=zeros(L,L);
        
%         fitnessgenerous=zeros(L,L);
%         rhogenerous=zeros(L,L);

        featen=zeros(L,L,N-1);
        geaten=zeros(L,L,N-1);

        fgettoeat=zeros(L,L,N-1);
        ggettoeat=zeros(L,L,N-1);
        
        fboth=zeros(L,L,N-1);
        gboth=zeros(L,L,N-1);
        
%         fgenerous=zeros(L,L,N-1);
%         ggenerous=zeros(L,L,N-1);
        
        parfor ind=1:(L*L*(N-1))
%         parfor ind=1:5
            [u,v,k]=ind2sub([L,L,N-1],ind);
            % rows are residents, columns are invaders, transpose to plot
            % normally
            resident=strats(u);
            invader=strats(v);

%             strategy=resident*ones(1,N);
%             strategy(N+1-(1:k))=invader;
%             [probeaten, probgettoeat, both, generous]=signalingevents_hypergeo_parallel(strategy,numsigs_permove,nummoves,m);

            toadd=zeros(m+1,1);
            for l=0:m
               toadd(l+1)=hygepdf(l,N-1,m,invader)*power(1-hygecdf(l,N-1,m,invader),k-1)*power(1-hygecdf(floor(l/invader*resident),N-1,m,resident),N-k);
            end
            probleast=sum(toadd);
            probeaten=zeros(1,N);
            probeaten(1:k)=probleast;
            probeaten((k+1):end)=(1-probleast*k)/(N-k);
            
            toadd=zeros(m+1,1);
            for l=0:m
%                toadd(l+1)=hygepdf(l,N-1,m,invader)*power(hygecdf(l-1,N-1,m,invader),k-1)*power(hygecdf(l/invader*resident,N-1,m,resident),N-k);
               toadd(l+1)=hygepdf(l,N-1,m,invader)*power(hygecdf(l-1,N-1,m,invader),k-1)*power(hygecdf(floor(l/invader*resident),N-1,m,resident),N-k);
            end
            probmost=sum(toadd);
            probgettoeat=zeros(1,N);
            probgettoeat(1:k)=probmost;
            probgettoeat((k+1):end)=(1-probmost*k)/(N-k);
            
            both=(1-probeaten).*probgettoeat;

            perfeaten=1-probeaten;
            featen(ind)=mean(perfeaten(N+1-(1:k)));
            geaten(ind)=mean(perfeaten(N+1-((k+1):N)));

            perfgettoeat=probgettoeat;
            fgettoeat(ind)=mean(perfgettoeat(N+1-(1:k)));
            ggettoeat(ind)=mean(perfgettoeat(N+1-((k+1):N)));
            
            perfboth=both;
            fboth(ind)=mean(perfboth(N+1-(1:k)));
            gboth(ind)=mean(perfboth(N+1-((k+1):N)));
            
%             perfgenerous=generous;
%             fgenerous(ind)=mean(perfgenerous(N+1-(1:k)));
%             ggenerous(ind)=mean(perfgenerous(N+1-((k+1):N)));

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
            
            fboth_now=fboth(where);
            gboth_now=gboth(where);
            
%             fgenerous_now=fgenerous(where);
%             ggenerous_now=ggenerous(where);

            fitnesseaten(ind)=featen_now(1)/geaten_now(1);
            fitnessgettoeat(ind)=fgettoeat_now(1)/ggettoeat_now(1);
            fitnessboth(ind)=fboth_now(1)/gboth_now(1);
%             fitnessgenerous(ind)=fgenerous_now(1)/ggenerous_now(1);

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
            
            ratio=gboth_now./fboth_now;
            tosum=ones(1,N-1);
            for k=1:(N-1)
                tosum(k)=prod(ratio(1:k));
            end
            rhoboth(ind)=1/(1+sum(tosum));
%             
%             ratio=ggenerous_now./fgenerous_now;
%             tosum=ones(1,N-1);
%             for k=1:(N-1)
%                 tosum(k)=prod(ratio(1:k));
%             end
%             rhogenerous(ind)=1/(1+sum(tosum));
        end  
        
        storefitnesseaten{q,p}=fitnesseaten;
        storerhoeaten{q,p}=rhoeaten;
        
        storefitnessgettoeat{q,p}=fitnessgettoeat;
        storerhogettoeat{q,p}=rhogettoeat;
        
        storefitnessboth{q,p}=fitnessboth;
        storerhoboth{q,p}=rhoboth;
        
%         storefitnessgenerous{q,p}=fitnessgenerous;
%         storerhogenerous{q,p}=rhogenerous;

        eqstratseaten=eq_strats(N,fitnesseaten,rhoeaten);
        ESSseaten{q,p}=eqstratseaten{2};

        eqstratsgettoeat=eq_strats(N,fitnessgettoeat,rhogettoeat);
        ESSsgettoeat{q,p}=eqstratsgettoeat{2};
        
        eqstratsboth=eq_strats(N,fitnessboth,rhoboth);
        ESSsboth{q,p}=eqstratsboth{2};
        
%         eqstratsgenerous=eq_strats(N,fitnessgenerous,rhogenerous);
%         ESSsgenerous{q,p}=eqstratsgenerous{2};
    end
end

t=toc(t0);
disp(t);

filename=strcat('/home/brush/schooling_consensus/ESS_4regimes_hypogeom_really','_nummoves=',num2str(nummoves),'_numpermove=',num2str(numsigs_permove),'.mat');
% save(filename,'strats','mvals','Tvals','storefitnesseaten','storerhoeaten','ESSseaten','storefitnessgettoeat','storerhogettoeat','ESSsgettoeat','storefitnessboth','storerhoboth','ESSsboth','storefitnessgenerous','storerhogenerous','ESSsgenerous');
save(filename,'strats','mvals','Tvals','storefitnesseaten','storerhoeaten','ESSseaten','storefitnessgettoeat','storerhogettoeat','ESSsgettoeat','storefitnessboth','storerhoboth','ESSsboth');

delete(dellapool);

exit ;
