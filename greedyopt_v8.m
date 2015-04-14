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

timesteps=100;
its=20;

evolution_eaten=zeros(its,N,timesteps);
H2s_eaten=zeros(its,timesteps);
H2sforced_eaten=zeros(its,timesteps);
corrlengths_eaten=zeros(its,timesteps);
corrlengthsforced_eaten=zeros(its,timesteps);
disconnected_eaten=zeros(its,timesteps);
rhos_eaten=zeros(its,timesteps);
rhosforced_eaten=zeros(its,timesteps);
changeorder_eaten=zeros(its,N,timesteps);

evolution_gettoeat=zeros(its,N,timesteps);
H2s_gettoeat=zeros(its,timesteps);
H2sforced_gettoeat=zeros(its,timesteps);
corrlengths_gettoeat=zeros(its,timesteps);
corrlengthsforced_gettoeat=zeros(its,timesteps);
disconnected_gettoeat=zeros(its,timesteps);
rhos_gettoeat=zeros(its,timesteps);
rhosforced_gettoeat=zeros(its,timesteps);
changeorder_gettoeat=zeros(its,N,timesteps);

evolution_both=zeros(its,N,timesteps);
H2s_both=zeros(its,timesteps);
H2sforced_both=zeros(its,timesteps);
corrlengths_both=zeros(its,timesteps);
corrlengthsforced_both=zeros(its,timesteps);
disconnected_both=zeros(its,timesteps);
rhos_both=zeros(its,timesteps);
rhosforced_both=zeros(its,timesteps);
changeorder_both=zeros(its,N,timesteps);

evolution_generous=zeros(its,N,timesteps);
H2s_generous=zeros(its,timesteps);
H2sforced_generous=zeros(its,timesteps);
corrlengths_generous=zeros(its,timesteps);
corrlengthsforced_generous=zeros(its,timesteps);
disconnected_generous=zeros(its,timesteps);
rhos_generous=zeros(its,timesteps);
rhosforced_generous=zeros(its,timesteps);
changeorder_generous=zeros(its,N,timesteps);

maxt_eaten=zeros(its,1);
maxt_gettoeat=zeros(its,1);
maxt_both=zeros(its,1);
maxt_generous=zeros(its,1);

% for ir=1:Nr
    radius=radvals(aid);
for num=1:its
    init=randi([1 N-1],1,N);
    t=1;
    timesteps2=timesteps;
    first_stop=0;
    while t<=timesteps2

    if t==1
        strategy=init;
    %     strategy=2*ones(1,N);
    %     strategy(1)=10; 
        evolution_eaten(num,:,t)=strategy;
        t=t+1;
    else
        strategy=evolution_eaten(num,:,t-1);
        strategy=reshape(strategy,N,1);

%         i=randsample(1:N,1);
%         for i=1:N
        changeorder=randsample(1:N,N,'false');
        changeorder_eaten(num,:,t)=changeorder;     
        for i=changeorder
            [probeaten, ~,~,~]=signalingevents_new_v2_parallel(strategy,numsigs_permove,nummoves,radius,b,T);
            perf=1-probeaten;
            
            choosestrat=zeros(1,3);
            choosestrat(2)=perf(i);

            if strategy(i)==N-1
                choosestrat(3)=perf(i);
            else 
                upone=strategy;
                upone(i)=strategy(i)+1;

                [probeaten, ~,~,~]=signalingevents_new_v2_parallel(upone,numsigs_permove,nummoves,radius,b,T);
                choosestrat(3)=1-probeaten(i);
            end
            
            if strategy(i)==0
                choosestrat(1)=perf(i);
            else
                downone=strategy;
                downone(i)=strategy(i)-1;

                [probeaten, ~,~,~]=signalingevents_new_v2_parallel(downone,numsigs_permove,nummoves,radius,b,T);
                choosestrat(1)=1-probeaten(i);   
            end 
            
            if choosestrat(1)>choosestrat(2)
                strategy(i)=strategy(i)-1;
            end
            if choosestrat(3)>choosestrat(2)
                strategy(i)=strategy(i)+1;
            end
            
        end

        evolution_eaten(num,:,t)=strategy;
        
        if sum(evolution_eaten(num,:,t)==evolution_eaten(num,:,t-1))==N
            evolution_eaten(num,:,(t+1):end)=repmat(col(strategy),1,timesteps-t);
            if first_stop==0
                t_eaten=t;
                first_stop=1;
                timesteps2=min(t_eaten+5,timesteps);
            end
            t=t+1;
        else
            t_eaten=t;
            t=t+1;
        end
%         if sum(evolution_eaten(num,:,t)==evolution_eaten(num,:,t-1))==N
%             evolution_eaten(num,:,(t+1):end)=repmat(col(strategy),1,timesteps-t);
%             t_eaten=t;
%             t=timesteps+1;
%         else
%             t_eaten=t;
%             t=t+1; 
%         end
    end
    end
    maxt_eaten(num)=t_eaten;

    t=1;
    timesteps2=timesteps;
    first_stop=0;
    while t<=timesteps2

    if t==1
        strategy=init;
    %     strategy=2*ones(1,N);
    %     strategy(1)=10; 
        evolution_gettoeat(num,:,t)=strategy;
        t=t+1;
    else
        strategy=evolution_gettoeat(num,:,t-1);
        strategy=reshape(strategy,N,1);

%         i=randsample(1:N,1);
%         for i=1:N
        changeorder=randsample(1:N,N,'false');
        changeorder_gettoeat(num,:,t)=changeorder;     
        for i=changeorder
            [~, probgettoeat,~,~]=signalingevents_new_v2_parallel(strategy,numsigs_permove,nummoves,radius,b,T);
            perf=probgettoeat;
            
            choosestrat=zeros(1,3);
            choosestrat(2)=perf(i);

            if strategy(i)==N-1
                choosestrat(3)=perf(i);
            else 
                upone=strategy;
                upone(i)=strategy(i)+1;

                [~, probgettoeat,~,~]=signalingevents_new_v2_parallel(upone,numsigs_permove,nummoves,radius,b,T);
                choosestrat(3)=probgettoeat(i);
            end

            if strategy(i)==0
                choosestrat(1)=perf(i);
            else
                downone=strategy;
                downone(i)=strategy(i)-1;

                [~, probgettoeat,~,~]=signalingevents_new_v2_parallel(downone,numsigs_permove,nummoves,radius,b,T);
                choosestrat(1)=1-probgettoeat(i);   
            end 

            if choosestrat(1)>choosestrat(2)
                strategy(i)=strategy(i)-1;
            end
            if choosestrat(3)>choosestrat(2)
                strategy(i)=strategy(i)+1;
            end

        end
        evolution_gettoeat(num,:,t)=strategy;
        if sum(evolution_gettoeat(num,:,t)==evolution_gettoeat(num,:,t-1))==N || max(evolution_gettoeat(num,:,t-1))==1
            evolution_gettoeat(num,:,(t+1):end)=repmat(col(strategy),1,timesteps-t);
            if first_stop==0
                t_gettoeat=t;
                first_stop=1;
                timesteps2=min(t_gettoeat+5,timesteps);
            end
            t=t+1;
        else 
            t_gettoeat=t;
            t=t+1;
        end
%         if sum(evolution_gettoeat(num,:,t)==evolution_gettoeat(num,:,t-1))==N
%             evolution_gettoeat(num,:,(t+1):end)=repmat(col(strategy),1,timesteps-t);
%             t_gettoeat=t;
%             t=timesteps+1;
%         else 
%             t_gettoeat=t;
%             t=t+1;
%         end
    end
    end
    maxt_gettoeat(num)=t_gettoeat;
    
    t=1;
    timesteps2=timesteps;
    first_stop=0;
    while t<=timesteps2

    if t==1
        strategy=init;
    %     strategy=2*ones(1,N);
    %     strategy(1)=10; 
        evolution_both(num,:,t)=strategy;
        t=t+1;
    else
        strategy=evolution_both(num,:,t-1);
        strategy=reshape(strategy,N,1);

%         i=randsample(1:N,1);
%         for i=1:N
        changeorder=randsample(1:N,N,'false');
        changeorder_both(num,:,t)=changeorder;     
        for i=changeorder
            [~,~,both, ~]=signalingevents_new_v2_parallel(strategy,numsigs_permove,nummoves,radius,b,T);
            perf=both;
            
            choosestrat=zeros(1,3);
            choosestrat(2)=perf(i);

            if strategy(i)==N-1
                choosestrat(3)=perf(i);
            else 
                upone=strategy;
                upone(i)=strategy(i)+1;

                [~,~,both, ~]=signalingevents_new_v2_parallel(upone,numsigs_permove,nummoves,radius,b,T);
                choosestrat(3)=both(i);
            end
            
            if strategy(i)==0
                choosestrat(1)=perf(i);
            else
                downone=strategy;
                downone(i)=strategy(i)-1;

                [~,~,both, ~]=signalingevents_new_v2_parallel(downone,numsigs_permove,nummoves,radius,b,T);
                choosestrat(1)=both(i);   
            end 
            
            if choosestrat(1)>choosestrat(2)
                strategy(i)=strategy(i)-1;
            end
            if choosestrat(3)>choosestrat(2)
                strategy(i)=strategy(i)+1;
            end
            
        end

        evolution_both(num,:,t)=strategy;
        
        if sum(evolution_both(num,:,t)==evolution_both(num,:,t-1))==N
            evolution_both(num,:,(t+1):end)=repmat(col(strategy),1,timesteps-t);
            if first_stop==0
                t_both=t;
                first_stop=1;
                timesteps2=min(t_both+5,timesteps);
            end
            t=t+1;
        else
            t_both=t;
            t=t+1;
        end
%         if sum(evolution_both(num,:,t)==evolution_both(num,:,t-1))==N
%             evolution_both(num,:,(t+1):end)=repmat(col(strategy),1,timesteps-t);
%             t_both=t;
%             t=timesteps+1;
%         else
%             t_both=t;
%             t=t+1; 
%         end
    end
    end
    maxt_both(num)=t_both;
    
    t=1;
    timesteps2=timesteps;
    first_stop=0;
    while t<=timesteps2

    if t==1
        strategy=init;
    %     strategy=2*ones(1,N);
    %     strategy(1)=10; 
        evolution_generous(num,:,t)=strategy;
        t=t+1;
    else
        strategy=evolution_generous(num,:,t-1);
        strategy=reshape(strategy,N,1);

%         i=randsample(1:N,1);
%         for i=1:N
        changeorder=randsample(1:N,N,'false');
        changeorder_generous(num,:,t)=changeorder;     
        for i=changeorder
            [~,~,~,generous]=signalingevents_new_v2_parallel(strategy,numsigs_permove,nummoves,radius,b,T);
            perf=generous;
            
            choosestrat=zeros(1,3);
            choosestrat(2)=perf(i);

            if strategy(i)==N-1
                choosestrat(3)=perf(i);
            else 
                upone=strategy;
                upone(i)=strategy(i)+1;

                [~,~,~,generous]=signalingevents_new_v2_parallel(upone,numsigs_permove,nummoves,radius,b,T);
                choosestrat(3)=generous(i);
            end
            
            if strategy(i)==0
                choosestrat(1)=perf(i);
            else
                downone=strategy;
                downone(i)=strategy(i)-1;

                [~,~,~,generous]=signalingevents_new_v2_parallel(downone,numsigs_permove,nummoves,radius,b,T);
                choosestrat(1)=generous(i);   
            end 
            
            if choosestrat(1)>choosestrat(2)
                strategy(i)=strategy(i)-1;
            end
            if choosestrat(3)>choosestrat(2)
                strategy(i)=strategy(i)+1;
            end
            
        end

        evolution_generous(num,:,t)=strategy;
        
        if sum(evolution_generous(num,:,t)==evolution_generous(num,:,t-1))==N
            evolution_generous(num,:,(t+1):end)=repmat(col(strategy),1,timesteps-t);
            if first_stop==0
                t_generous=t;
                first_stop=1;
                timesteps2=min(t_generous+5,timesteps);
            end
            t=t+1;
        else
            t_generous=t;
            t=t+1;
        end
%         if sum(evolution_generous(num,:,t)==evolution_generous(num,:,t-1))==N
%             evolution_generous(num,:,(t+1):end)=repmat(col(strategy),1,timesteps-t);
%             t_generous=t;
%             t=timesteps+1;
%         else
%             t_generous=t;
%             t=t+1; 
%         end
    end
    end
    maxt_generous(num)=t_generous;
end

% end
clear num t init timesteps2 strategy perf probeaten probgettoeat both generous
delete(dellapool);
clear dellapool dellacluster

filename=strcat('/home/brush/schooling_consensus/greedyopt','_T=',num2str(T),'_nummoves=',num2str(nummoves),'_numpermove=',num2str(numsigs_permove),'_rad=',num2str(radius),'_timesteps=',num2str(timesteps),'.mat');
% save(filename,'evolution_eaten','evolution_gettoeat','changeorder_eaten','changeorder_gettoeat','disconnected_gettoeat','radvals','maxt_eaten','maxt_gettoeat')
save(filename)
'Greedy opt is done. Group props are not.' %#ok<NOPTS>

% for ir=1:Nr
for num=1:its
    for t=1:maxt_eaten(num)
        strategy_eaten=evolution_eaten(num,:,t);
        [meanH2, meanH2_forced, corrlength, corrlength_forced, disconnectedcount, meanrho, meanrho_forced]=groupprops(strategy_eaten,numsigs_permove,nummoves,radius,b,T);
        H2s_eaten(num,t)=meanH2;
        H2sforced_eaten(num,t)=meanH2_forced;
        corrlengths_eaten(num,t)=corrlength;
        corrlengthsforced_eaten(num,t)=corrlength_forced;
        rhos_eaten(num,t)=meanrho;
        rhosforced_eaten(num,t)=meanrho_forced;
        disconnected_eaten(num,t)=disconnectedcount;
    end
    H2s_eaten(num,(t+1):end)=meanH2;
    H2sforced_eaten(num,(t+1):end)=meanH2_forced;
    corrlengths_eaten(num,(t+1):end)=corrlength;
    corrlengthsforced_eaten(num,(t+1):end)=corrlength_forced;
    disconnected_eaten(num,(t+1):end)=disconnectedcount;
    rhos_eaten(num,(t+1):end)=meanrho;
    rhosforced_eaten(num,(t+1):end)=meanrho_forced;
    for t=1:maxt_gettoeat(num)
        strategy_gettoeat=evolution_gettoeat(num,:,t);
        [meanH2, meanH2_forced, corrlength, corrlength_forced, disconnectedcount, meanrho, meanrho_forced]=groupprops(strategy_gettoeat,numsigs_permove,nummoves,radius,b,T);
        H2s_gettoeat(num,t)=meanH2;
        H2sforced_gettoeat(num,t)=meanH2_forced;
        corrlengths_gettoeat(num,t)=corrlength;
        corrlengthsforced_gettoeat(num,t)=corrlength_forced;
        disconnected_gettoeat(num,t)=disconnectedcount;
        rhos_gettoeat(num,t)=meanrho;
        rhosforced_gettoeat(num,t)=meanrho_forced;
    end
    H2s_gettoeat(num,(t+1):end)=meanH2;
    H2sforced_gettoeat(num,(t+1):end)=meanH2_forced;
    corrlengths_gettoeat(num,(t+1):end)=corrlength;
    corrlengthsforced_gettoeat(num,(t+1):end)=corrlength_forced;
    disconnected_gettoeat(num,(t+1):end)=disconnectedcount;
    rhos_gettoeat(num,(t+1):end)=meanrho;
    rhosforced_gettoeat(num,(t+1):end)=meanrho_forced;
    for t=1:maxt_both(num)
        strategy_both=evolution_both(num,:,t);
        [meanH2, meanH2_forced, corrlength, corrlength_forced, disconnectedcount, meanrho, meanrho_forced]=groupprops(strategy_both,numsigs_permove,nummoves,radius,b,T);
        H2s_both(num,t)=meanH2;
        H2sforced_both(num,t)=meanH2_forced;
        corrlengths_both(num,t)=corrlength;
        corrlengthsforced_both(num,t)=corrlength_forced;
        disconnected_both(num,t)=disconnectedcount;
        rhos_both(num,t)=meanrho;
        rhosforced_both(num,t)=meanrho_forced;
    end
    H2s_both(num,(t+1):end)=meanH2;
    H2sforced_both(num,(t+1):end)=meanH2_forced;
    corrlengths_both(num,(t+1):end)=corrlength;
    corrlengthsforced_both(num,(t+1):end)=corrlength_forced;
    disconnected_both(num,(t+1):end)=disconnectedcount;
    rhos_both(num,(t+1):end)=meanrho;
    rhosforced_both(num,(t+1):end)=meanrho_forced;
    for t=1:maxt_generous(num)
        strategy_generous=evolution_generous(num,:,t);
        [meanH2, meanH2_forced, corrlength, corrlength_forced, disconnectedcount, meanrho, meanrho_forced]=groupprops(strategy_generous,numsigs_permove,nummoves,radius,b,T);
        H2s_generous(num,t)=meanH2;
        H2sforced_generous(num,t)=meanH2_forced;
        corrlengths_generous(num,t)=corrlength;
        corrlengthsforced_generous(num,t)=corrlength_forced;
        disconnected_generous(num,t)=disconnectedcount;
        rhos_generous(num,t)=meanrho;
        rhosforced_generous(num,t)=meanrho_forced;
    end
    H2s_generous(num,(t+1):end)=meanH2;
    H2sforced_generous(num,(t+1):end)=meanH2_forced;
    corrlengths_generous(num,(t+1):end)=corrlength;
    corrlengthsforced_generous(num,(t+1):end)=corrlength_forced;
    disconnected_generous(num,(t+1):end)=disconnectedcount;
    rhos_generous(num,(t+1):end)=meanrho;
    rhosforced_generous(num,(t+1):end)=meanrho_forced;
end
% end

clear num t meanH2 meanH2_forced corrlength corrlength_forced disconnectedcount meanrho meanrho_forced

filename=strcat('/home/brush/schooling_consensus/greedyopt','_T=',num2str(T),'_nummoves=',num2str(nummoves),'_numpermove=',num2str(numsigs_permove),'_rad=',num2str(radius),'_timesteps=',num2str(timesteps),'.mat');
% save(filename,'evolution_eaten','H2s_eaten','H2sforced_eaten','corrlengths_eaten','corrlengthsforced_eaten','rhos_eaten','rhosforced_eaten','evolution_gettoeat','H2s_gettoeat','H2sforced_gettoeat','corrlengths_gettoeat','corrlengthsforced_gettoeat','rhos_gettoeat','rhosforced_gettoeat','evolution_both','H2s_both','H2sforced_both','corrlengths_both','corrlengthsforced_both','rhos_both','rhosforced_both','evolution_generous','H2s_generous','H2sforced_generous','corrlengths_generous','corrlengthsforced_generous','rhos_generous','rhosforced_generous','disconnected_eaten','disconnected_gettoeat','disconnected_both','disconnected_generous','radvals','maxt_eaten','maxt_gettoeat','maxt_both','maxt_generous','changeorder_eaten','changeorder_gettoeat','changeorder_both','changeorder_generous','radius','b','T','nummoves','numsigs_permove','N')
save(filename)
'Greedy opt is done. Group props are too.' %#ok<NOPTS>



exit ;
