timesteps=3;

evolution_eaten=zeros(N,(timesteps-1)*N+1);
groupconsensus_eaten=zeros(1,(timesteps-1)*N+1);
groupconsensusforced_eaten=zeros(1,(timesteps-1)*N+1);
corrlengths_eaten=zeros(1,(timesteps-1)*N+1);
corrlengthsforced_eaten=zeros(1,(timesteps-1)*N+1);
disconnected_eaten=zeros(1,(timesteps-1)*N+1);
changeorder_eaten=zeros(N,timesteps);

evolution_gettoeat=zeros(N,(timesteps-1)*N+1);
groupconsensus_gettoeat=zeros(1,(timesteps-1)*N+1);
groupconsensusforced_gettoeat=zeros(1,(timesteps-1)*N+1);
corrlengths_gettoeat=zeros(1,(timesteps-1)*N+1);
corrlengthsforced_gettoeat=zeros(1,(timesteps-1)*N+1);
disconnected_gettoeat=zeros(1,(timesteps-1)*N+1);
changeorder_gettoeat=zeros(N,timesteps);


k=15;
% init=[14*ones(k,1);17*ones(N-k,1)];
init=randi([1 N-1],1,N);


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


% t=1;
% while t<=timesteps
% 
% if t==1
%     strategy=init;
% %     strategy=2*ones(1,N);
% %     strategy(1)=10; 
%     evolution_eaten(:,t)=strategy;
%     t=t+1;
% else
%     strategy=evolution_eaten(:,(t-2)*N+1);
%     strategy=reshape(strategy,N,1);
% 
% %         [probeaten, ~]=signalingevents_parallel(strategy,numsigs_permove,nummoves,radius,b,T);
% %         perf=1-probeaten;
% %         newstrategy=strategy;
% 
% %         i=randsample(1:N,1);
% %         for i=1:N
%     changeorder=randsample(1:N,N,'false');
%     for count=1:N
%         i=changeorder(count);
%         [probeaten, ~]=signalingevents_parallel(strategy,numsigs_permove,nummoves,radius,b,T);
%         perf=1-probeaten;
% 
%         choosestrat=zeros(1,3);
%         choosestrat(2)=perf(i);
% 
%         if strategy(i)==N-1
%             choosestrat(3)=perf(i);
%         else 
%             upone=strategy;
%             upone(i)=strategy(i)+1;
% 
%             [probeaten, ~]=signalingevents_parallel(upone,numsigs_permove,nummoves,radius,b,T);
%             choosestrat(3)=1-probeaten(i);
%         end
% 
%         if strategy(i)==1
%             choosestrat(1)=perf(i);
%         else
%             downone=strategy;
%             downone(i)=strategy(i)-1;
% 
%             [probeaten, ~]=signalingevents_parallel(downone,numsigs_permove,nummoves,radius,b,T);
%             choosestrat(1)=1-probeaten(i);   
%         end 
% 
%         if choosestrat(1)>choosestrat(2)
%             strategy(i)=strategy(i)-1;
%         end
%         if choosestrat(3)>choosestrat(2)
%             strategy(i)=strategy(i)+1;
%         end
%         evolution_eaten(:,(t-2)*N+count+1)=strategy;
%     end
%     changeorder_eaten(:,t)=changeorder;
%     if sum(evolution_eaten(:,(t-1)*N+1)==evolution_eaten(:,(t-2)*N+1))==N
%         evolution_eaten(:,((t-1)*N+2):end)=repmat(col(strategy),1,(timesteps-t)*N);
%         t_eaten=t;
%         t=timesteps+1;
%     else
%         t_eaten=t;
%         t=t+1;
%     end
% end
% 
% end

t=1;
while t<=timesteps

if t==1
    strategy=init;
%     strategy=2*ones(1,N);
%     strategy(1)=10; 
    evolution_gettoeat(:,t)=strategy;
    t=t+1;
else
    strategy=evolution_gettoeat(:,(t-2)*N+1);
    strategy=reshape(strategy,N,1);

%         [~, probgettoeat]=signalingevents_parallel(strategy,numsigs_permove,nummoves,radius,b,T);
%         perf=probgettoeat;
%         newstrategy=strategy;

%         i=randsample(1:N,1);
%         for i=1:N
    changeorder=randsample(1:N,N,'false');
%     for i=changeorder
    for count=1:N
        i=changeorder(count);
        [~, probgettoeat]=signalingevents_parallel(strategy,numsigs_permove,nummoves,radius,b,T);
        perf=probgettoeat;

        choosestrat=zeros(1,3);
        choosestrat(2)=perf(i);

        if strategy(i)==N-1
            choosestrat(3)=perf(i);
        else 
            upone=strategy;
            upone(i)=strategy(i)+1;

            [~, probgettoeat]=signalingevents_parallel(upone,numsigs_permove,nummoves,radius,b,T);
            choosestrat(3)=probgettoeat(i);
        end

        if strategy(i)==1
            choosestrat(1)=perf(i);
        else
            downone=strategy;
            downone(i)=strategy(i)-1;

            [~, probgettoeat]=signalingevents_parallel(downone,numsigs_permove,nummoves,radius,b,T);
            choosestrat(1)=1-probgettoeat(i);   
        end 

        if choosestrat(1)>choosestrat(2)
            strategy(i)=strategy(i)-1;
        end
        if choosestrat(3)>choosestrat(2)
            strategy(i)=strategy(i)+1;
        end
        evolution_gettoeat(:,(t-2)*N+count+1)=strategy;
    end
    changeorder_gettoeat(:,t)=changeorder;     
%     if sum(evolution_gettoeat(:,(t-1)*N+1)==evolution_gettoeat(:,(t-2)*N+1))==N
%         evolution_gettoeat(:,((t-1)*N+2):end)=repmat(col(strategy),1,(timesteps-t)*N);
%         t_gettoeat=t;
%         t=timesteps+1;
%     else 
%         t_gettoeat=t;
        t=t+1;
%     end
end
end

%%
i=1;
strategy=evolution_eaten(i,:,end);
strategy=sort(strategy);
[probeaten,~]=signalingevents(strategy,numsigs_permove,nummoves,radius,b,T);
