
groupconsensus_eaten=zeros(its,timesteps);
groupconsensusforced_eaten=zeros(its,timesteps);
corrlengths_eaten=zeros(its,timesteps);
corrlengthsforced_eaten=zeros(its,timesteps);
disconnected_eaten=zeros(its,timesteps);

groupconsensus_gettoeat=zeros(its,timesteps);
groupconsensusforced_gettoeat=zeros(its,timesteps);
corrlengths_gettoeat=zeros(its,timesteps);
corrlengthsforced_gettoeat=zeros(its,timesteps);
disconnected_gettoeat=zeros(its,timesteps);

for num=1:its
    for t=1:maxt_eaten(num)
        strategy_eaten=evolution_eaten(num,:,t);
        [meanH2, meanH2_forced, corrlength, corrlength_forced, disconnectedcount]=groupprops(strategy_eaten,numsigs_permove,nummoves,radius,b,T);
        groupconsensus_eaten(num,t)=meanH2;
        groupconsensusforced_eaten(num,t)=meanH2_forced;
        corrlengths_eaten(num,t)=corrlength;
        corrlengthsforced_eaten(num,t)=corrlength_forced;
        disconnected_eaten(num,t)=disconnectedcount;
    end
    groupconsensus_eaten(num,(t+1):end)=meanH2;
    groupconsensusforced_eaten(num,(t+1):end)=meanH2_forced;
    corrlengths_eaten(num,(t+1):end)=corrlength;
    corrlengthsforced_eaten(num,(t+1):end)=corrlength_forced;
    disconnected_eaten(num,(t+1):end)=disconnectedcount;
    for t=1:maxt_gettoeat(num)
        strategy_gettoeat=evolution_gettoeat(num,:,t);
        [meanH2, meanH2_forced, corrlength, corrlength_forced, disconnectedcount]=groupprops(strategy_gettoeat,numsigs_permove,nummoves,radius,b,T);
        groupconsensus_gettoeat(num,t)=meanH2;
        groupconsensusforced_gettoeat(num,t)=meanH2_forced;
        corrlengths_gettoeat(num,t)=corrlength;
        corrlengthsforced_gettoeat(num,t)=corrlength_forced;
        disconnected_gettoeat(num,t)=disconnectedcount;
    end
    groupconsensus_gettoeat(num,(t+1):end)=meanH2;
    groupconsensusforced_gettoeat(num,(t+1):end)=meanH2_forced;
    corrlengths_gettoeat(num,(t+1):end)=corrlength;
    corrlengthsforced_gettoeat(num,(t+1):end)=corrlength_forced;
    disconnected_gettoeat(num,(t+1):end)=disconnectedcount;
end