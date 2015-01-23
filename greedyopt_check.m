num=10;
f=find(col(sum(abs(diff(evolution_gettoeat(1,:,:),[],3)),2))==0,1,'first');
strategy=evolution_gettoeat(num,:,f);
[~, perf]=signalingevents_parallel(strategy,numsigs_permove,nummoves,radius,b,T);
newstrategy=strategy;
for i=1:N
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
        newstrategy(i)=strategy(i)-1;
    end
    if choosestrat(3)>choosestrat(2)
        newstrategy(i)=strategy(i)+1;
    end

end
