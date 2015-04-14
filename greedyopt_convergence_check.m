evolution_eaten=zeros(its,N,timesteps);
changeorder_eaten=zeros(its,N,timesteps);
num=1;
timesteps=10;

% init=randi([1 N-1],1,N);
init=19*ones(N,1);
init(6)=11;init(11)=16;

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
        [probeaten, ~,~,~]=signalingevents_new_v2(strategy,numsigs_permove,nummoves,radius,b,T);
        perf=1-probeaten;

        choosestrat=zeros(1,3);
        choosestrat(2)=perf(i);

        if strategy(i)==N-1
            choosestrat(3)=perf(i);
        else 
            upone=strategy;
            upone(i)=strategy(i)+1;

            [probeaten, ~,~,~]=signalingevents_new_v2(upone,numsigs_permove,nummoves,radius,b,T);
            choosestrat(3)=1-probeaten(i);
        end

        if strategy(i)==0
            choosestrat(1)=perf(i);
        else
            downone=strategy;
            downone(i)=strategy(i)-1;

            [probeaten, ~,~,~]=signalingevents_new_v2(downone,numsigs_permove,nummoves,radius,b,T);
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