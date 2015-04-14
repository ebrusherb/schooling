function [costs_to_low,costs_to_high] = costs_by_category(N,strats,invaderstrats,nummoves,radius,b,T)

mainstrats = strats;
Ls = length(mainstrats);

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
    strategy = strat*ones(1,N);
    strategy(invader) = invaderstrat;

    scores = zeros(N,nummoves);

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

        keepMtot(:,:,i,k)=M;
        keepreceivertot(i,k)=receiver;
        keepreceiverstot(:,i,k)=allreceivers;
        
    end
    
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
        minscorer(rows(look),i) = 1/size(look,1);
    %     minscorer(rows(cols==i),i) = 1/numsigs_tot;
    end
    
    keepminscorestot(:,:,k)=minscorer;
        

end



cats=cell(4,Ls);
cat_indicators=[0 0;0 1;1 0;1 1];
onerec = cell(Ls,1);
morerec = cell(Ls,1);
numonerec = zeros(Ls,1);
nummorerec = zeros(Ls,1);
minscores_bycat = cell(5,Ls);

for k=1:Ls
    numrec = sum(keepreceiverstot(:,:,k),1);
    resrec = find(keepreceiverstot(invader-1,:,k)==1);
    invrec = find(keepreceiverstot(invader,:,k)==1);
    onerec{k} = setdiff(find(numrec==1),union(resrec,invrec));
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
    
%     morerec{k} = setdiff(find(numrec>1),union(resrec,invrec));
%     nummorerec(k) = length(morerec{k});
    
    cats{5,k}=setdiff(find(numrec>1),union(resrec,invrec));
    nummorerec(k)=length(cats{5,k});
    minscores_bycat{5,k}=keepminscorestot([v(low) v(high)],cats{5,k},k);
     
end

costs_to_low=zeros(5,Ls);
costs_to_high=zeros(5,Ls);
numcats=zeros(5,Ls);
meanbycat=zeros(2,5,Ls);

for k=1:Ls    
    numcats_now=zeros(2,5);
    for i=1:5
    
        numcats(i,k)=length(cats{i,k});
        numcats_now(:,i)=length(cats{i,k});
        if ~isempty(cats{i,k})
            meanbycat(:,i,k)=mean(minscores_bycat{i,k},2);
        end
        
    end
    todivideby=(numonerec(k)+nummorerec(k))*ones(2,5);
    costs_now=meanbycat(:,:,k).*numcats_now./todivideby;
    costs_to_low(:,k)=costs_now(1,:);
    costs_to_high(:,k)=costs_now(2,:);
end
