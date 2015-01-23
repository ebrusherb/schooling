% numworkers=str2num(getenv('PROCS'));
% dellacluster=parcluster('local');
% dellacluster.JobStorageLocation=strcat('/scratch/network/brush/tmp/',getenv('SLURM_JOB_ID'));
% dellapool=parpool(dellacluster, numworkers) ;

numsigs_permove=1;
nummoves=10;

N=10;
b=1;
radius=.1;
T=1;

timesteps=5;
its=2;

evolution_eaten=zeros(its,N,timesteps);
groupspeed_eaten=zeros(its,timesteps);
groupconsensus_eaten=zeros(its,timesteps);
corrlengths_eaten=zeros(its,timesteps);
disconnected_eaten=zeros(its,timesteps);

evolution_gettoeat=zeros(its,N,timesteps);
groupspeed_gettoeat=zeros(its,timesteps);
groupconsensus_gettoeat=zeros(its,timesteps);
corrlengths_gettoeat=zeros(its,timesteps);
disconnected_gettoeat=zeros(its,timesteps);

max_eaten=zeros(1,its);
max_gettoeat=zeros(1,its);

%strategy evolution
for num=1:its
    init=randi([1 N-1],1,N);

'predation regime'
t=1;
while t<=timesteps
% for t=1:timesteps    
    
if t==1
    strategy=init;
%     strategy=2*ones(1,N);
%     strategy(1)=10; 
    evolution_eaten(num,:,t)=strategy;
    t=t+1;
else
    strategy=evolution_eaten(num,:,t-1);
    strategy=reshape(strategy,1,N);

    [probeaten, ~]=signalingevents_parallel(strategy,numsigs_permove,nummoves,radius,b,T);
    perf=1-probeaten;
    newstrategy=strategy;

    for i=1:N
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

    end
        evolution_eaten(num,:,t)=newstrategy;
%     if sum(col(strategy)==col(newstrategy))<N
%         evolution_eaten(num,:,t)=newstrategy;
%         t
        max_eaten(num)=t;
        t=t+1
%     else 
        evolution_eaten(num,:,t:end)=repmat(col(newstrategy),1,timesteps-t+1);
%         t=timesteps+1
%         timesteps+1
%     end
end
end

'resource regime'
% t=1;
% while t<=timesteps
%     
% if t==1
%     strategy=init;
% %     strategy=2*ones(1,N);
% %     strategy(1)=10; 
%     evolution_gettoeat(num,:,t)=strategy;
% else
%     strategy=evolution_gettoeat(num,:,t-1);
%     strategy=reshape(strategy,1,N);
% 
%     [~, probgettoeat]=signalingevents_parallel(strategy,numsigs_permove,nummoves,radius,b,T);
%     perf=probgettoeat;
%     newstrategy=strategy;
% 
%     for i=1:N
%         choosestrat=zeros(1,3);
%         choosestrat(2)=perf(i);
% 
%         if strategy(i)==N-1
%             choosestrat(3)=perf(i);
%         else 
%             upone=strategy;
%             upone(i)=strategy(i)+1;
% 
%             [~, probgettoeat]=signalingevents_parallel(upone,numsigs_permove,nummoves,radius,b,T);
%             choosestrat(3)=probgettoeat(i);
%         end
% 
%         if strategy(i)==1
%             choosestrat(1)=perf(i);
%         else
%             downone=strategy;
%             downone(i)=strategy(i)-1;
% 
%             [~, probgettoeat]=signalingevents_parallel(downone,numsigs_permove,nummoves,radius,b,T);
%             choosestrat(1)=1-probgettoeat(i);   
%         end 
% 
%         if choosestrat(1)>choosestrat(2)
%             newstrategy(i)=strategy(i)-1;
%         end
%         if choosestrat(3)>choosestrat(2)
%             newstrategy(i)=strategy(i)+1;
%         end
% 
%     end
%     if sum(col(strategy)==col(newstrategy))<N
%         evolution_gettoeat(num,:,t)=newstrategy;
%         max_gettoeat(num)=t
%         t=t+1;        
%     else
%         evolution_gettoeat(num,:,t:end)=repmat(col(newstrategy),1,timesteps-t+1);
%         t=timesteps+1
%     end
% end
% end
end

% filename=strcat('/home/brush/schooling_consensus/greedyopt','_rad=',num2str(radius),'_T=',num2str(T),'_nummoves=',num2str(nummoves),'_numpermove=',num2str(numsigs_permove),'timesteps=',num2str(timesteps),'.mat');
% save(filename,'evolution_eaten','groupspeed_eaten','groupconsensus_eaten','corrlengths_eaten','evolution_gettoeat','groupspeed_gettoeat','groupconsensus_gettoeat','corrlengths_gettoeat','disconnected_eaten','disconnected_gettoeat')
% 
% 'Greedy opt is done. Group props are not.' %#ok<NOPTS>
% 
% %group properties
% for num=1:its
%     for t=1:max_eaten(num)
%         strategy_eaten=evolution_eaten(num,:,t);
%         [meanlambda, meanH2, meancorrlength, disconnectedcount]=groupprops_parallel(strategy_eaten,numsigs_permove,nummoves,radius,b,T);
%         groupspeed_eaten(num,t)=meanlambda;
%         groupconsensus_eaten(num,t)=meanH2;
%         corrlengths_eaten(num,t)=meancorrlength;
%         disconnected_eaten(num,t)=disconnectedcount;
%     end
%         groupspeed_eaten(num,(max_eaten(num)+1):end)=repmat(meanlambda,timesteps-max_eaten(num),1);
%         groupconsensus_eaten(num,(max_eaten(num)+1):end)=repmat(meanH2,timesteps-max_eaten(num),1);
%         corrlengths_eaten(num,(max_eaten(num)+1):end)=repmat(meancorrlength,timesteps-max_eaten(num),1);
%         disconnected_eaten(num,(max_eaten(num)+1):end)=repmat(disconnectedcount,timesteps-max_eaten(num),1);
%     for t=1:max_gettoeat(num)
%         strategy_gettoeat=evolution_gettoeat(num,:,t);
%         [meanlambda, meanH2, meancorrlength, disconnectedcount]=groupprops_parallel(strategy_gettoeat,numsigs_permove,nummoves,radius,b,T);
%         groupspeed_gettoeat(num,t)=meanlambda;
%         groupconsensus_gettoeat(num,t)=meanH2;
%         corrlengths_gettoeat(num,t)=meancorrlength;
%         disconnected_gettoeat(num,t)=disconnectedcount;
%     end
%         groupspeed_gettoeat(num,(max_gettoeat(num)+1):end)=repmat(meanlambda,timesteps-max_gettoeat(num),1);
%         groupconsensus_gettoeat(num,(max_gettoeat(num)+1):end)=repmat(meanH2,timesteps-max_gettoeat(num),1);
%         corrlengths_gettoeat(num,(max_gettoeat(num)+1):end)=repmat(meancorrlength,timesteps-max_gettoeat(num),1);
%         disconnected_gettoeat(num,(max_gettoeat(num)+1):end)=repmat(disconnectedcount,timesteps-max_gettoeat(num),1);
% end
% 
% filename=strcat('/home/brush/schooling_consensus/greedyopt','_rad=',num2str(radius),'_T=',num2str(T),'_nummoves=',num2str(nummoves),'_numpermove=',num2str(numsigs_permove),'timesteps=',num2str(timesteps),'.mat');
% save(filename,'evolution_eaten','groupspeed_eaten','groupconsensus_eaten','corrlengths_eaten','evolution_gettoeat','groupspeed_gettoeat','groupconsensus_gettoeat','corrlengths_gettoeat','disconnected_eaten','disconnected_gettoeat')
% 
% 'Greedy opt is done. Group props are too.' %#ok<NOPTS>
% 
% delete(dellapool);
% 
% exit ;
