%%

N=20;
its=10000;
temp=2;

randits=10;
finalstrats=zeros(N,randits);
H2vovertime=zeros(randits,its+1);
corrlengthovertime=zeros(randits,its+1);

numsigs_permove=1;
nummoves=100;
radius=.1;
b=1;
T=1;

parfor r=1:randits
    H2vec=zeros(its+1,1);
    corrvec=zeros(its+1,1);
    strategy=zeros(N,its+1);
    perf=zeros(N,its+1);
    mutant=zeros(1,its+1);

    recordofchanges=zeros(1,its+1);

    strategy(:,1)=2*randi([1 floor((N-1)/2)],N,1);
    [meanH2, meancorrlength, probeaten, ~]=signalingevents(strategy(:,1),numsigs_permove,nummoves,radius,b,T);
    perf(:,1)=probeaten;
    H2vec(1)=meanH2;
    corrvec(1)=meancorrlength;
    for time=(1:its)+1
        strategy(:,time)=strategy(:,time-1);

        changer=randi([1 N],1);
        mutant(time)=changer;
        test=[reshape(strategy(:,time-1),N,1) reshape(strategy(:,time-1),N,1)];
        v=exp(-abs((2:2:floor((N-1)/2))-test(changer,1))/temp);
        probs=v/sum(v);
        test(changer,2)=randchoice(2:2:floor((N-1)/2),probs);
        [meanH2, meancorrlength, probeaten, ~]=signalingevents(test(:,2),numsigs_permove,nummoves,radius,b,T);
        changed=probeaten;
        performance=[perf(:,time-1) changed];
        [~,w]=min(performance(changer,:));
        strategy(changer,time)=test(changer,w);
        H2vec(time)=meanH2;
        corrvec(time)=meancorrlength;
        if w==2 && strategy(changer,time-1)~=strategy(changer,time)
            recordofchanges(time)=1;
        end

        [~, ~, probeaten, ~]=signalingevents(strategy(:,time),numsigs_permove,nummoves,radius,b,T);
        perf(:,time)=probeaten;
    end

    finalstrats(:,r)=strategy(:,end);
    H2overtime(r,:)=H2vec;
    corrlengthovertime(r,:)=corrvec;
end
