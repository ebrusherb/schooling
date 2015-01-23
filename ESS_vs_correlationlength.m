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

%% find prob eaten for number of informed birds and number of informed neighbors
N=20;
radius = 1.5;
nummoves=1000;
T = 1;
strats = [16 18];
invaderstrats = fliplr(strats);

maxrec = zeros(1,2);
invader_counts = zeros(N,N,2);
invader_minscores = cell(N,N,2);
resident_counts = zeros(N,N,2);
resident_minscores = cell(N,N,2);
keepscorestot = zeros(N,nummoves,2);
keepMtot = zeros(N,N,nummoves,2);
keepdtot = zeros(N,N,nummoves,2);
keepreceivertot = zeros(nummoves,2);
keepreceiverstot = zeros(N,nummoves,2);
keepminscorestot = zeros(N,nummoves,2);

probeaten = zeros(N,2);

for k = 1:2
    strat = strats(k);
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


        keepscorestot(:,i,k)=v;
        keepMtot(:,:,i,k)=M;
        keepreceivertot(i,k)=receiver;
        keepreceiverstot(:,i,k)=allreceivers;
        
    end
    
    keepM=reshape(keepM,N,N,nummoves);
    keepreceivers = reshape(keepreceivers,[],nummoves);
    numrec = sum(keepreceivers,1);
    maxrec(k) = max(numrec);
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

    payingattention = keepM;
    payingattention(keepM == -1) = 1;
    tosignal = payingattention;
    tonaive = payingattention;
    for i = 1:nummoves
        tosignal(:,:,i) = payingattention(:,:,i).*repmat(col(keepreceivers(:,i))',N,1);
        tonaive(:,:,i) = payingattention(:,:,i).*repmat(col(~keepreceivers(:,i))',N,1);
    end
    sumtosignal = reshape(sum(tosignal,2),N,nummoves);
    sumtonaive = reshape(sum(tonaive,2),N,nummoves);

    normsumtosignal = sigfig(sumtosignal.*repmat(strategy',1,nummoves));
    
    normsumtosignal_res = normsumtosignal(resident,:);
    minscorer_res = minscorer(resident,:);
    keepreceivers_res = keepreceivers(resident,:);
        
    for i = 1:N
        for j = 1:N
        finv = find(numrec==i);
        finv2 = find(normsumtosignal(invader,:)==j-1);
        finv3 = find(keepreceivers(invader,:)==1);
        finv4 = setdiff(intersect(finv,finv2),finv3);
        invader_counts(i,j,k) = length(finv4);
        invader_minscores{i,j,k} = col(minscorer(invader,finv4));
        
        fres = find(repmat(numrec,N-1,1)==i);
        fres2 = find(normsumtosignal_res==j-1);
        fres3 = find(keepreceivers_res==1);
        fres4 = setdiff(intersect(fres,fres2),fres3);
        resident_counts(i,j,k) = length(fres4);
        resident_minscores{i,j,k} = col(minscorer_res(fres4));
        end
    end 

end

%%
totalprob_invader = zeros(N,N,2);
totalprob_resident = zeros(N,N,2);

for i = 1:N
    for j = 1:N
        for k = 1:2
            if invader_counts(i,j,k)>0
                totalprob_invader(i,j,k) = sum(invader_minscores{i,j,k})/(1/nummoves)/invader_counts(i,j,k);
            end
        end
    end
end
for i = 1:N
    for j = 1:N
        for k = 1:2
            if resident_counts(i,j,k)>0
                totalprob_resident(i,j,k) = sum(resident_minscores{i,j,k})/(1/nummoves)/resident_counts(i,j,k);
            end
        end
    end
end

maxrec = max(maxrec);
%%
cats=cell(4,2);

cat_indicators=[0 0;0 1;1 0;1 1];
onerec = cell(2,1);
numonerec = zeros(2,1);
minscores_bycat = cell(4,2);

for k=1:2
    numrec = sum(keepreceiverstot(:,:,k),1);
    onerec{k} = find(numrec==1);
    numonerec(k) = length(onerec{k});
    
    for i=onerec{k}
        receiver = keepreceivertot(i,k);
        for j=1:4
        [~,low]=min([strats(k) invaderstrats(k)]);
        [~,high]=max([strats(k) invaderstrats(k)]);
        v=[invader-1 invader];
        if receiver~=invader-1 && receiver~=invader && (keepMtot(v(high),receiver,i,k)>0)==cat_indicators(j,1) && (keepMtot(v(low),receiver,i,k)>0)==cat_indicators(j,2)
            cats{j,k}=[cats{j,k},i];
            minscores_bycat{j,k}=horzcat(minscores_bycat{j,k},keepminscorestot([v(low) v(high)],i,k));
        end
        end
    end
end

meanbycat=zeros(2,4,2);
varbycat=zeros(2,4,2);
numcats=zeros(2,4,2);

for i=1:4
    for k=1:2
        numcats(1,i,k)=length(cats{i,k});
        numcats(2,i,k)=length(cats{i,k});
        if numcats(1,i,k)>0
            meanbycat(:,i,k)=mean(minscores_bycat{i,k},2)*nummoves;
            varbycat(:,i,k)=var(minscores_bycat{i,k},[],2);
        end
    end
end

todivideby=zeros(2,4,2);
todivideby(:,:,1)=numonerec(1);
todivideby(:,:,2)=numonerec(2);

probs=sum(meanbycat.*numcats./todivideby,2);
%%
mycolormap=cbrewer('seq','Reds',10);

subplot(2,4,1)
imagesc(1:maxrec,0:(maxrec),transpose(totalprob_resident(1:maxrec,1:(maxrec+1),1)))


set(gca,'ydir','normal')
set(gca,'xtick',1:maxrec,'xticklabel',1:maxrec,'ytick',0:maxrec,'yticklabel',0:maxrec,'FontName',fontname,'FontSize',labfontsz)
colormap(mycolormap)
caxis manual
caxis([0 1])
ylab=ylabel('Informed neighbors','FontName',fontname,'FontSize',textfontsz);
xlabel('Informed birds','FontName',fontname,'FontSize',textfontsz)
% title(strcat('Resident = ',num2str(strats(1)),', Invader = ',num2str(invaderstrats(1))),'FontName',fontname,'FontSize',textfontsz);
title(strcat('Res Prob: ','Resident = ',num2str(strats(1)),', Invader = ',num2str(invaderstrats(1))),'FontName',fontname,'FontSize',textfontsz);
set(gca,'FontName',fontname,'FontSize',labfontsz)

subplot(2,4,2)
imagesc(1:maxrec,0:(maxrec),transpose(totalprob_invader(1:maxrec,1:(maxrec+1),1)))
v=get(gca,'Position');
set(gca,'ydir','normal')
set(gca,'xtick',1:maxrec,'xticklabel',1:maxrec,'ytick',0:maxrec,'yticklabel',0:maxrec,'FontName',fontname,'FontSize',labfontsz)
colormap(mycolormap)
caxis manual
caxis([0 1])
ylab=ylabel('Informed neighbors','FontName',fontname,'FontSize',textfontsz);
xlabel('Informed birds','FontName',fontname,'FontSize',textfontsz)
% title(strcat('Resident = ',num2str(strats(1)),', Invader = ',num2str(invaderstrats(1))),'FontName',fontname,'FontSize',textfontsz);
title(strcat('Invader Prob: ','Resident = ',num2str(strats(1)),', Invader = ',num2str(invaderstrats(1))),'FontName',fontname,'FontSize',textfontsz);
set(gca,'FontName',fontname,'FontSize',labfontsz)

subplot(2,4,3)
imagesc(1:maxrec,0:(maxrec),transpose(totalprob_resident(1:maxrec,1:(maxrec+1),1).*resident_counts(1:maxrec,1:(maxrec+1),1)/(N-1)))
% imagesc(1:maxrec,0:(maxrec),transpose(totalprob_invader(1:maxrec,1:(maxrec+1),1)./totalprob_resident(1:maxrec,1:(maxrec+1),1)))
% imagesc(1:maxrec,0:(maxrec),transpose(totalprob_invader(1:maxrec,1:(maxrec+1),1).*invader_counts(1:maxrec,1:(maxrec+1),1)./(totalprob_resident(1:maxrec,1:(maxrec+1),1).*resident_counts(1:maxrec,1:(maxrec+1),1))))
colorbar
set(gca,'ydir','normal')
set(gca,'xtick',1:maxrec,'xticklabel',1:maxrec,'ytick',0:maxrec,'yticklabel',0:maxrec,'FontName',fontname,'FontSize',labfontsz)
colormap(mycolormap)
caxis manual
% caxis([0 1])

subplot(2,4,4)
imagesc(1:maxrec,0:(maxrec),transpose(totalprob_invader(1:maxrec,1:(maxrec+1),1).*invader_counts(1:maxrec,1:(maxrec+1),1)))
% imagesc(1:maxrec,0:(maxrec),transpose(totalprob_invader(1:maxrec,1:(maxrec+1),1)./totalprob_resident(1:maxrec,1:(maxrec+1),1)))
% imagesc(1:maxrec,0:(maxrec),transpose(totalprob_invader(1:maxrec,1:(maxrec+1),1).*invader_counts(1:maxrec,1:(maxrec+1),1)./(totalprob_resident(1:maxrec,1:(maxrec+1),1).*resident_counts(1:maxrec,1:(maxrec+1),1))))
colorbar
set(gca,'ydir','normal')
set(gca,'xtick',1:maxrec,'xticklabel',1:maxrec,'ytick',0:maxrec,'yticklabel',0:maxrec,'FontName',fontname,'FontSize',labfontsz)
colormap(mycolormap)
caxis manual
% caxis([0 1])

subplot(2,4,5)
imagesc(1:maxrec,0:maxrec,transpose(totalprob_resident(1:maxrec,1:(maxrec+1),2)))
v=get(gca,'Position');
set(gca,'ydir','normal')
set(gca,'xtick',1:maxrec,'xticklabel',1:maxrec,'ytick',0:maxrec,'yticklabel',0:maxrec,'FontName',fontname,'FontSize',labfontsz)
colormap(mycolormap)
caxis manual
caxis([0 1])
ylab=ylabel('Informed neighbors','FontName',fontname,'FontSize',textfontsz);
xlabel('Informed birds','FontName',fontname,'FontSize',textfontsz)
% title(strcat('Resident = ',num2str(strats(1)),', Invader = ',num2str(invaderstrats(1))),'FontName',fontname,'FontSize',textfontsz);
title(strcat('Res Prob: ','Resident = ',num2str(strats(2)),', Invader = ',num2str(invaderstrats(2))),'FontName',fontname,'FontSize',textfontsz);
set(gca,'FontName',fontname,'FontSize',labfontsz)

subplot(2,4,6)
imagesc(1:maxrec,0:maxrec,transpose(totalprob_invader(1:maxrec,1:(maxrec+1),2)))
v=get(gca,'Position');
set(gca,'ydir','normal')
set(gca,'xtick',1:maxrec,'xticklabel',1:maxrec,'ytick',1:maxrec,'yticklabel',1:maxrec,'FontName',fontname,'FontSize',labfontsz)
colormap(mycolormap)
caxis manual
caxis([0 1])
ylab=ylabel('Informed neighbors','FontName',fontname,'FontSize',textfontsz);
xlabel('Informed birds','FontName',fontname,'FontSize',textfontsz)
% title(strcat('Resident = ',num2str(strats(1)),', Invader = ',num2str(invaderstrats(1))),'FontName',fontname,'FontSize',textfontsz);
title(strcat('Invader Prob: ','Resident = ',num2str(strats(2)),', Invader = ',num2str(invaderstrats(2))),'FontName',fontname,'FontSize',textfontsz);
set(gca,'FontName',fontname,'FontSize',labfontsz)

subplot(2,4,7)
imagesc(1:maxrec,0:(maxrec),transpose(totalprob_resident(1:maxrec,1:(maxrec+1),2).*(resident_counts(1:maxrec,1:(maxrec+1),2)/(N-1))))
% imagesc(1:maxrec,0:(maxrec),transpose(totalprob_invader(1:maxrec,1:(maxrec+1),2)./(totalprob_resident(1:maxrec,1:(maxrec+1),2))))
% imagesc(1:maxrec,0:(maxrec),transpose(totalprob_invader(1:maxrec,1:(maxrec+1),2).*invader_counts(1:maxrec,1:(maxrec+1),2)./(totalprob_resident(1:maxrec,1:(maxrec+1),2).*resident_counts(1:maxrec,1:(maxrec+1),2))))
set(gca,'ydir','normal')
set(gca,'xtick',1:maxrec,'xticklabel',1:maxrec,'ytick',0:maxrec,'yticklabel',0:maxrec,'FontName',fontname,'FontSize',labfontsz)
colormap(mycolormap)
caxis manual
colorbar
% caxis([0 1])

subplot(2,4,8)
imagesc(1:maxrec,0:(maxrec),transpose(totalprob_invader(1:maxrec,1:(maxrec+1),2).*(invader_counts(1:maxrec,1:(maxrec+1),2))))
% imagesc(1:maxrec,0:(maxrec),transpose(totalprob_invader(1:maxrec,1:(maxrec+1),2)./(totalprob_resident(1:maxrec,1:(maxrec+1),2))))
% imagesc(1:maxrec,0:(maxrec),transpose(totalprob_invader(1:maxrec,1:(maxrec+1),2).*invader_counts(1:maxrec,1:(maxrec+1),2)./(totalprob_resident(1:maxrec,1:(maxrec+1),2).*resident_counts(1:maxrec,1:(maxrec+1),2))))
set(gca,'ydir','normal')
set(gca,'xtick',1:maxrec,'xticklabel',1:maxrec,'ytick',0:maxrec,'yticklabel',0:maxrec,'FontName',fontname,'FontSize',labfontsz)
colormap(mycolormap)
caxis manual
colorbar
% caxis([0 1])
%%
close all
k=1;
f1=find(keepminscorertot(invader,:,k)*nummoves>.95);
i=i+1;
M=reshape(keepMtot(:,:,f1(i),k),N,N);
d=reshape(keepdtot(:,:,f1(i),k),N,N);
receiver=keepreceivertot(f1(i),k);
beta = zeros(N,1);
allreceivers = d(receiver,:)<= radius;
beta(allreceivers) = b;
v = real(expected_spin(M,T,beta));
numsubset=keepreceiverstot(:,f1,k);
numtot=sum(numsubset,1);
[distvec,corrvec,corrlength]=correlationlength_mat_single(M,d,b,radius,receiver);
corrs=paircorrelations(M,beta);

hold on
% plot(d(receiver,~allreceivers),v(~allreceivers),'o','LineWidth',2)
% plot(d(receiver,invader),v(invader),'or','LineWidth',2)
plot(corrs(receiver,~allreceivers),v(~allreceivers),'o','LineWidth',2)
plot(corrs(receiver,invader),v(invader),'or','LineWidth',2)
% plot(d(receiver,~allreceivers),corrs(receiver,~allreceivers),'o','LineWidth',2)
% plot(d(receiver,invader),corrs(receiver,invader),'or','LineWidth',2)
set(gcf,'Position',[850 386 560 420])

%%
strats = [16 17];
invaderstrats = fliplr(strats);
strategies=zeros(N,2);
strategies(:,1)=strats(1);
strategies(10,1)=invaderstrats(1);
strategies(:,2)=strats(2);
strategies(9,2)=invaderstrats(2);

testits=1000;
problowleftbehind=zeros(N,testits,2);
scoreslowleftbehind=zeros(N,testits,2);
for i=1:testits
    for k=1:2
        strategy=strategies(:,k);
        positions = unifrnd(0,1,N,2);
        d = squareform(pdist(positions));

        keepdtot(:,:,i,k) = d;

        M = zeros(N);
        for ind = 1:N
            [~, order] = sort(d(ind,:));
            neighbors = order(2:strategy(ind)+1);
            M(ind,neighbors) = 1/strategy(ind);
        end
        M(1:N+1:end) = -1; 

        receiver=find(M(9,:)==0,1,'first');

        beta = zeros(N,1);
        beta(receiver) = b;
        v = real(expected_spin(M,T,beta));
        scoreslowleftbehind(:,i,k)=v;
        f = find(v==min(v));
        probeaten = zeros(N,1);
        probeaten(f) = 1/length(f)/testits;
        problowleftbehind(:,i,k)=probeaten;
    end
end

probhighdiluted=zeros(N,testits,2);
scoreshighdiluted=zeros(N,testits,2);

for i=1:testits
    for k=1:2
        strategy=strategies(:,k);
        positions = unifrnd(0,1,N,2);
        d = squareform(pdist(positions));

        keepdtot(:,:,i,k) = d;

        M = zeros(N);
        for ind = 1:N
            [~, order] = sort(d(ind,:));
            neighbors = order(2:strategy(ind)+1);
            M(ind,neighbors) = 1/strategy(ind);
        end
        M(1:N+1:end) = -1; 

        receiver=find(M(10,:)>0,1,'first');

        beta = zeros(N,1);
        beta(receiver) = b;
        v = real(expected_spin(M,T,beta));
        scoreshighdiluted(:,i,k)=v;
        f = find(v==min(v));
        probeaten = zeros(N,1);
        probeaten(f) = 1/length(f)/testits;
        probhighdiluted(:,i,k)=probeaten;
    end
end
%%
f1=find(probhighdiluted(10,:,1)*testits==1);
f2=find(probhighdiluted(10,:,2)*testits==1);


