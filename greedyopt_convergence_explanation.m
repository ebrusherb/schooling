timesteps=100;

evolution_eaten=zeros(N,timesteps);
groupconsensus_eaten=zeros(timesteps);
groupconsensusforced_eaten=zeros(timesteps);
corrlengths_eaten=zeros(timesteps);
corrlengthsforced_eaten=zeros(timesteps);
disconnected_eaten=zeros(timesteps);


k=15;
init=[14*ones(k,1);17*ones(N-k,1)];
% [probeaten, ~]=signalingevents_parallel(strategy,numsigs_permove,nummoves,radius,b,T);
% [mean(1-probeaten(1:k)) mean(1-probeaten(k+1:end))]
% 
% strategy_upone=strategy;
% strategy_upone(k+1)=strategy(end)+1;
% 
% strategy_downone=strategy;
% strategy_downone(k+1)=strategy(end)-1;
% 
% [probeaten_upone, ~]=signalingevents_parallel(strategy_upone,numsigs_permove,nummoves,radius,b,T);
% [probeaten_downone, ~]=signalingevents_parallel(strategy_downone,numsigs_permove,nummoves,radius,b,T);

t=1;
    while t<=timesteps

    if t==1
        strategy=init;
    %     strategy=2*ones(1,N);
    %     strategy(1)=10; 
        evolution_eaten(:,t)=strategy;
        t=t+1;
    else
        strategy=evolution_eaten(:,t-1);
        strategy=reshape(strategy,1,N);

        [probeaten, ~]=signalingevents_parallel(strategy,numsigs_permove,nummoves,radius,b,T);
        perf=1-probeaten;
        newstrategy=strategy;
        
        i=randsample(1:N,1);
%         for i=1:N
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

%         end
        evolution_eaten(:,t)=newstrategy;
        if sum(evolution_eaten(:,t)==evolution_eaten(:,t-1))==N
            evolution_eaten(:,(t+1):end)=repmat(col(newstrategy),1,timesteps-t);
%             t_eaten=t;
            t=timesteps+1;
        else
%             t_eaten=t;
            t=t+1; 
        end
    end
    end