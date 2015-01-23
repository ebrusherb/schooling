%% load ESS data

ESSeaten_total=cell(2,1);
ESSgettoeat_total=cell(2,1);

load /Users/eleanorbrush/Desktop/ESS_nummoves=1000_numpermove=1_fullstorage.mat

ESSeaten_total{1}=ESSseaten;
ESSgettoeat_total{1}=ESSsgettoeat;

load /Users/eleanorbrush/Desktop/ESS_nummoves=10_numpermove=100.mat

ESSeaten_total{2}=ESSseaten;
ESSgettoeat_total{2}=ESSsgettoeat;

load /Users/eleanorbrush/Desktop/greedyopt_rad=0.1_T=1_nummoves=10_numpermove=100.mat

load /Users/eleanorbrush/Desktop/homogengroupprops.mat
%% THIS IS PRETTY SCREWED UP: I COPIED AND PASTED AND COPIED AGAIN AND NOWS ITS ALL A MESS.
N=20;
radius = .1;
b = 1;
nummoves=1000;
T = 1;
% strats=18*ones(17,1);
mainstrats = 3:17;
% strats=[14:17];
Ls = length(mainstrats);
% invaderstrats = 1:17;
% invaderstrats = 18*ones(Ls,1);
invaderstrats = mainstrats+1;
% invaderstrats = fliplr(strats);


keepMtot = zeros(N,N,nummoves,Ls);
keepdtot = zeros(N,N,nummoves,Ls);
keepreceivertot = zeros(nummoves,Ls);
keepreceiverstot = zeros(N,nummoves,Ls);
keepminscorestot = zeros(N,nummoves,Ls);

probeaten = zeros(N,2);

for k = 1:Ls
    strat = mainstrats(k);
    invader = 10;
    resident = [1:(invader-1) (invader+1):N];
    invaderstrat = invaderstrats(k);
    strategy = [strat*ones(1,N)];
    strategy(invader) = invaderstrat;

    scores = zeros(N,nummoves);
    keepM = zeros(N,N,nummoves);
    keepreceivers = zeros(N,nummoves);

    % parfor i=1:nummoves
    for i = 1:nummoves
        positions = unifrnd(0,1,N,2);
        d = squareform(pdist(positions));
        
        keepdtot(:,:,i,k) = d;

        M = zeros(N);
        for ind = 1:N
            [~, order] = sort(d(ind,:));
            neighbors = order(2:strategy(ind)+1);
            M(ind,neighbors) = 1/strategy(ind);
        end
        M(1:N+1:end) = -1; %sets diagonal equal to -1
    %     M = M(perm,perm);
        
        beta = zeros(N,1);
        receiver = randsample(N,1,'true');
        allreceivers = d(receiver,:)<= radius;
        beta(allreceivers) = b;
        v = real(expected_spin(M,T,beta));

        scores(:,i) = v;
        keepM(:,:,i)= M;
        keepreceivers(:,i) = (allreceivers);


        
        keepMtot(:,:,i,k)=M;
        keepreceivertot(i,k)=receiver;
        keepreceiverstot(:,i,k)=allreceivers;
        
    end
    
    keepM=reshape(keepM,N,N,nummoves);
    keepreceivers = reshape(keepreceivers,[],nummoves);
    numrec = sum(keepreceivers,1);
    
    scores = reshape(scores,N,[]);
    scores = sigfig(scores,4);
    [c,o] = sort(scores,1);

    % minscorer = o(1,:);

    [minvals,~] = min(scores);
    minmat = repmat(minvals,N,1);
    [rows,cols] = find(abs(scores-minmat)<0.00001);
    minscorer = zeros(N,size(scores,2));
    
    for i = 1:size(scores,2)
        look = find(cols==i);
        minscorer(rows(look),i) = 1/size(look,1)/nummoves;
    %     minscorer(rows(cols==i),i) = 1/numsigs_tot;
    end
    
    keepminscorestot(:,:,k)=minscorer;
    probeaten(:,k) = sum(minscorer,2);     

end

cats=cell(4,Ls);

cat_indicators=[0 0;0 1;1 0;1 1];
onerec = cell(Ls,1);
numonerec = zeros(Ls,1);
minscores_bycat = cell(4,Ls);

for k=1:Ls
    numrec = sum(keepreceiverstot(:,:,k),1);
    onerec{k} = find(numrec==1);
    numonerec(k) = length(onerec{k});
    
    for i=onerec{k}
        receiver = keepreceivertot(i,k);
        for j=1:4
        [~,low]=min([mainstrats(k) invaderstrats(k)]);
        [~,high]=max([mainstrats(k) invaderstrats(k)]);
        v=[invader-1 invader];
        if receiver~=invader-1 && receiver~=invader && (keepMtot(v(high),receiver,i,k)>0)==cat_indicators(j,1) && (keepMtot(v(low),receiver,i,k)>0)==cat_indicators(j,2)
            cats{j,k}=[cats{j,k},i];
            minscores_bycat{j,k}=horzcat(minscores_bycat{j,k},keepminscorestot([v(low) v(high)],i,k));
        end
        end
    end
end

costs_to_high=zeros(4,Ls);
numcats=zeros(4,Ls);
meanbycat=zeros(2,4,Ls);

for k=1:Ls    
    numcats_now=zeros(2,4);
    for i=1:4
    
        numcats(i,k)=length(cats{i,k});
        numcats_now(:,i)=length(cats{i,k});
        if ~isempty(cats{i,k})
            meanbycat(:,i,k)=mean(minscores_bycat{i,k},2)*nummoves;
        end
        todivideby=numonerec(k)*ones(2,4);
        costs=meanbycat(:,:,k).*numcats_now./todivideby;
        costs_to_high(:,k)=costs(2,:);
    end
end

% strats=18*ones(17,1);
mainstrats = 3:17;
% strats=[14:17];
Ls = length(mainstrats);
% invaderstrats = 1:17;
% invaderstrats = 2*ones(Ls,1);
invaderstrats = mainstrats-1;
% invaderstrats = fliplr(strats);


keepMtot = zeros(N,N,nummoves,Ls);
keepdtot = zeros(N,N,nummoves,Ls);
keepreceivertot = zeros(nummoves,Ls);
keepreceiverstot = zeros(N,nummoves,Ls);
keepminscorestot = zeros(N,nummoves,Ls);

probeaten = zeros(N,2);

for k = 1:Ls
    strat = mainstrats(k);
    invader = 10;
    resident = [1:(invader-1) (invader+1):N];
    invaderstrat = invaderstrats(k);
    strategy = [strat*ones(1,N)];
    strategy(invader) = invaderstrat;

    scores = zeros(N,nummoves);
    keepM = zeros(N,N,nummoves);
    keepreceivers = zeros(N,nummoves);

    % parfor i=1:nummoves
    for i = 1:nummoves
        positions = unifrnd(0,1,N,2);
        d = squareform(pdist(positions));
        
        keepdtot(:,:,i,k) = d;

        M = zeros(N);
        for ind = 1:N
            [~, order] = sort(d(ind,:));
            neighbors = order(2:strategy(ind)+1);
            M(ind,neighbors) = 1/strategy(ind);
        end
        M(1:N+1:end) = -1; %sets diagonal equal to -1
    %     M = M(perm,perm);
        
        beta = zeros(N,1);
        receiver = randsample(N,1,'true');
        allreceivers = d(receiver,:)<= radius;
        beta(allreceivers) = b;
        v = real(expected_spin(M,T,beta));

        scores(:,i) = v;
        keepM(:,:,i)= M;
        keepreceivers(:,i) = (allreceivers);

        keepMtot(:,:,i,k)=M;
        keepreceivertot(i,k)=receiver;
        keepreceiverstot(:,i,k)=allreceivers;
        
    end
    
    keepM=reshape(keepM,N,N,nummoves);
    keepreceivers = reshape(keepreceivers,[],nummoves);
    numrec = sum(keepreceivers,1);
   
    scores = reshape(scores,N,[]);
    scores = sigfig(scores,4);
    [c,o] = sort(scores,1);

    % minscorer = o(1,:);

    [minvals,~] = min(scores);
    minmat = repmat(minvals,N,1);
    [rows,cols] = find(abs(scores-minmat)<0.00001);
    minscorer = zeros(N,size(scores,2));
    
    for i = 1:size(scores,2)
        look = find(cols==i);
        minscorer(rows(look),i) = 1/size(look,1)/nummoves;
    %     minscorer(rows(cols==i),i) = 1/numsigs_tot;
    end
    
    keepminscorestot(:,:,k)=minscorer;
    probeaten(:,k) = sum(minscorer,2);     

end

cats=cell(4,Ls);

cat_indicators=[0 0;0 1;1 0;1 1];
onerec = cell(Ls,1);
numonerec = zeros(Ls,1);
minscores_bycat = cell(4,Ls);

for k=1:Ls
    numrec = sum(keepreceiverstot(:,:,k),1);
    resrec = find(keepreceivertot(:,k)==invader-1);
    invrec = find(keepreceivertot(:,k)==invader);
    onerec{k} = setdiff(find(numrec==1),intersect(resrec,invrec));
    numonerec(k) = length(onerec{k});
    
    for i=onerec{k}
        receiver = keepreceivertot(i,k);
        for j=1:4
        [~,low]=min([mainstrats(k) invaderstrats(k)]);
        [~,high]=max([mainstrats(k) invaderstrats(k)]);
        v=[invader-1 invader];
        if receiver~=invader-1 && receiver~=invader && (keepMtot(v(high),receiver,i,k)>0)==cat_indicators(j,1) && (keepMtot(v(low),receiver,i,k)>0)==cat_indicators(j,2)
            cats{j,k}=[cats{j,k},i];
            minscores_bycat{j,k}=horzcat(minscores_bycat{j,k},keepminscorestot([v(low) v(high)],i,k));
        end
        end
    end
end

costs_to_low=zeros(4,Ls);
numcats=zeros(4,Ls);
meanbycat=zeros(2,4,Ls);

for k=1:Ls    
    numcats_now=zeros(2,4);
    for i=1:4
    
        numcats(i,k)=length(cats{i,k});
        numcats_now(:,i)=length(cats{i,k});
        if ~isempty(cats{i,k})
            meanbycat(:,i,k)=mean(minscores_bycat{i,k},2)*nummoves;
        end
        todivideby=numonerec(k)*ones(2,4);
        costs=meanbycat(:,:,k).*numcats_now./todivideby;
        costs_to_low(:,k)=costs(1,:);
    end
end

%%
figure
mycols1=cbrewer('seq','YlOrRd',Ls);
mycols2=cbrewer('seq','YlGnBu',Ls);
stratlabs1=cell(1,Ls);
for i=1:Ls
    stratlabs1{i}=num2str(i);
end
stratlabs2=cell(1,Ls);
for i=1:Ls
    stratlabs2{i}=num2str(i+2);
end
plots1=zeros(1,Ls);
plots2=zeros(1,Ls);

subplot(1,2,1)
hold on
for i=1:Ls
    plots1(i)=plot(costs_to_high(:,i),'-o','Color',mycols1(i,:),'LineWidth',2); 
end
hold off
% set(gca,'ylim',[0 .6])
leg=legend(plots1,stratlabs1);
legend('boxoff');
set(leg,'position',[.2 .5 .04 .3])

subplot(1,2,2)
hold on
for i=1:Ls
    plots2(i)=plot(costs_to_low(:,i),'-o','Color',mycols2(i,:),'LineWidth',2);
end
hold off
% set(gca,'ylim',[0 .6])
leg=legend(plots2,stratlabs2);
legend('boxoff');
set(leg,'position',[.6 .5 .04 .3])

%%
figure
hold on

for i=1:Ls
%     plot(numcats(:,i),'-o','Color',mycols1(i,:),'LineWidth',2)
    plot(1:4,meanbycat(2,:,i),'-o','Color',mycols1(i,:),'LineWidth',2)
end
