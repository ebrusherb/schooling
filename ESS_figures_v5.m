N = 20;
nummoves = 1000;
numsigs_permove = 1;
b = 1;

fontname='Times';
textfontsz=12;
labfontsz=9;
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

%% find prob eaten for number of informed birds and number of informed neighbors
radius = .1;
T = 1;
strats = [16 17];
invaderstrats = fliplr(strats);

maxrec = zeros(1,2);
invader_counts = zeros(N,N,2);
invader_minscores = cell(N,N,2);
resident_counts = zeros(N,N,2);
resident_minscores = cell(N,N,2);

probeaten = zeros(N,2);

for k = 1:2
    strat = strats(k);
    invader = 10;
    resident = [1:(invader-1) (invader+1):N];
    invaderstrat = invaderstrats(k);
    strategy = [strat*ones(1,N)];
    strategy(invader) = invaderstrat;

    numsigs_tot = numsigs_permove*nummoves;

    scores = zeros(N,numsigs_permove,nummoves);
    keepM = zeros(N,N,numsigs_permove,nummoves);
    keepreceivers = zeros(N,numsigs_permove,nummoves);

    % parfor i=1:nummoves
    for i = 1:nummoves
        positions = unifrnd(0,1,N,2);
        d = squareform(pdist(positions));

        M = zeros(N);
        for ind = 1:N
            [~, order] = sort(d(ind,:));
            neighbors = order(2:strategy(ind)+1);
            M(ind,neighbors) = 1/strategy(ind);
        end
        M(1:N+1:end) = -1; %sets diagonal equal to -1
    %     M = M(perm,perm);
        receivers = randsample(N,numsigs_permove,'true');
        for j = 1:numsigs_permove
            beta = zeros(N,1);
            receiver = receivers(j);
            allreceivers = d(receiver,:)<= radius;
            beta(allreceivers) = b;
            v = real(expected_spin(M,T,beta));

            scores(:,j,i) = v;
            keepM(:,:,j,i)= M;
            keepreceivers(:,j,i) = (allreceivers);

        end

    end
    
    keepM=reshape(keepM,20,20,nummoves);
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
        minscorer(rows(look),i) = 1/size(look,1)/numsigs_tot;
    %     minscorer(rows(cols==i),i) = 1/numsigs_tot;
    end
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
        
    for i = 1:20
        for j = 1:N
        finv = find(numrec==i);
        finv2 = find(normsumtosignal(invader,:)==j);
        finv3 = find(keepreceivers(invader,:)==1);
        finv4 = setdiff(intersect(finv,finv2),finv3);
        invader_counts(i,j,k) = length(finv4);
        invader_minscores{i,j,k} = col(minscorer(invader,finv4));
        
        fres = find(repmat(numrec,N-1,1)==i);
        fres2 = find(normsumtosignal_res==j);
        fres3 = find(keepreceivers_res==1);
        fres4 = setdiff(intersect(fres,fres2),fres3);
        resident_counts(i,j,k) = length(fres4);
        resident_minscores{i,j,k} = col(minscorer_res(fres4));
        end
    end 

end


totalprob_invader = zeros(N,N,2);
totalprob_resident = zeros(N,N,2);

for i = 1:N
    for j = 1:N
        for k = 1:2
            if invader_counts(i,j,k)>0
                totalprob_invader(i,j,k) = sum(invader_minscores{i,j,k})/.001/invader_counts(i,j,k);
            end
        end
    end
end
for i = 1:N
    for j = 1:N
        for k = 1:2
            if resident_counts(i,j,k)>0
                totalprob_resident(i,j,k) = sum(resident_minscores{i,j,k})/.001/resident_counts(i,j,k);
            end
        end
    end
end

maxrec = max(maxrec);
%% find number of receivers for different radius
radstohist = [.1 .5 1.2];
T = 1;
Nr = length(radstohist);
numrecmat = zeros(Nr,nummoves);

for u = 1:Nr
    radius = radstohist(u);

    numsigs_tot = numsigs_permove*nummoves;
    keepreceivers = zeros(N,numsigs_permove,nummoves);
    
    % parfor i=1:nummoves
    for i = 1:nummoves
        positions = unifrnd(0,1,N,2);
        d = squareform(pdist(positions));

        M = zeros(N);
        for ind = 1:N
            [~, order] = sort(d(ind,:));
            neighbors = order(2:strategy(ind)+1);
            M(ind,neighbors) = 1/strategy(ind);
        end
        M(1:N+1:end) = -1; %sets diagonal equal to -1
    %     M = M(perm,perm);
        receivers = randsample(N,numsigs_permove,'true');
        for j = 1:numsigs_permove
            beta = zeros(N,1);
            receiver = receivers(j);
            allreceivers = d(receiver,:)<= radius;
            
            keepreceivers(:,j,i) = (allreceivers);

        end

    end
    
    numrecmat(u,:) = sum(keepreceivers,1);
end
%% ESS strategies and prob eaten comparison together
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=.75*w;
set(gcf,'Units','inches');
set(gcf,'Position',[11.5 3 w h]);
fontname='TimesNewRoman';
textfontsz=12;
labfontsz=10;
markersz=3;
lw=2;
axislw=1;

n=size(ESSseaten,2);
k=1;
t=3;

mycolormap = cbrewer('seq','Reds',11);

ESSseaten=ESSeaten_total{k};
ESSsgettoeat=ESSgettoeat_total{k};

subplot(3,2,1:2)
hold on

xdiv=[.25 1.05];
set(gca,'xlim',[-.05 1.45]);
xlim=get(gca,'xlim');
plot(xdiv(1)*ones(2,1),ylim,'--k','LineWidth',1)
plot(xdiv(2)*ones(2,1),ylim,'--k','LineWidth',1)
fill([xlim(1) xdiv(1) xdiv(1) xlim(1)],[0 0 N N],mycolormap(2,:),'EdgeColor','none')
fill([xdiv(1) xdiv(2) xdiv(2) xdiv(1)],[0 0 N N],mycolormap(5,:),'EdgeColor','none')
fill([xdiv(2) xlim(2) xlim(2) xdiv(2)],[0 0 N N],mycolormap(8,:),'EdgeColor','none')
for i=1:n
    l=size(ESSseaten{t,i},2);
    p1=plot(radvals(i)*ones(l,1),ESSseaten{t,i},'o','Color','k','MarkerSize',markersz,'LineWidth',lw);
    l=size(ESSsgettoeat{t,i},2);
    if l>0
        p2=plot(radvals(i)*ones(l,1),ESSsgettoeat{t,i},'s','Color','k','MarkerSize',markersz,'LineWidth',lw);
    end
end
xlab=xlabel('Radius of signal','FontName',fontname,'FontSize',textfontsz);
ylab=ylabel('ESS','FontName',fontname,'FontSize',textfontsz);
v=get(gca,'Position');
set(gca,'Position',[.1 v(2) .67 v(4)]) 
v=get(xlab,'Position');
set(xlab,'Position',[v(1) -4.5 v(3)]);
v=get(ylab,'Position');
set(ylab,'Position',[xlim(1)-.06*diff(xlim) v(2:3)]);
ylim=get(gca,'Ylim');
plot(xlim(1)*ones(2,1),ylim,'k','LineWidth',axislw);
plot(xlim,ylim(1)*ones(2,1),'k','LineWidth',axislw);

leg=legend([p1 p2],{'Predation','Resources'},'FontName',fontname,'FontSize',labfontsz);
legend('boxoff')
v=get(leg,'Position');
set(leg,'Position',[.77 .77 v(3) v(4)])

set(gca,'FontName',fontname,'FontSize',labfontsz)

[x,y]=data2norm(-.05,ylim(2));
annotation('textbox',[x-.09 1.03*y .01 .04],'String','(a)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

subplot(3,2,5)
imagesc(1:maxrec,1:maxrec,transpose(totalprob_invader(1:maxrec,1:maxrec,1)))
v=get(gca,'Position');
set(gca,'Position',[.1 .1 .32 .2])
set(gca,'ydir','normal')
set(gca,'xtick',1:maxrec,'xticklabel',1:maxrec,'ytick',1:maxrec,'yticklabel',1:maxrec,'FontName',fontname,'FontSize',labfontsz)
colormap(mycolormap)
caxis manual
caxis([0 1])
ylab=ylabel('Informed neighbors','FontName',fontname,'FontSize',textfontsz);
xlabel('Informed birds','FontName',fontname,'FontSize',textfontsz)
% title(strcat('Resident = ',num2str(strats(1)),', Invader = ',num2str(invaderstrats(1))),'FontName',fontname,'FontSize',textfontsz);
v=get(ylab,'Position');
set(ylab,'Position',[v(1)+.1 v(2)-.5 v(3)])

ylim=get(gca,'Ylim');
[x,y]=data2norm(0.5,ylim(2));
annotation('textbox',[x-.09 1.05*y .01 .04],'String','(c)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

subplot(3,2,6)
imagesc(1:maxrec,1:maxrec,transpose(totalprob_invader(1:maxrec,1:maxrec,2)))
v=get(gca,'Position');
set(gca,'Position',[.55 .1 .32 .2])
set(gca,'ydir','normal')
set(gca,'xtick',1:maxrec,'xticklabel',1:maxrec,'ytick',1:maxrec,'yticklabel',1:maxrec,'FontName',fontname,'FontSize',labfontsz)
colormap(mycolormap)
caxis manual
caxis([0 1])
c=colorbar;
v=get(c,'Position');
set(c,'Position',[.9 v(2)+.01 .02 .2],'Fontname',fontname,'FontSize',labfontsz)
% ylim=get(c,'YLabel');
% set(ylim,'String','Probability of being eaten','Rotation',270,'Position',[9 .5 1.001],'FontName',fontname,'FontSize',textfontsz);
ylab=ylabel('Informed neighbors','FontName',fontname,'FontSize',textfontsz);
xlabel('Informed birds','FontName',fontname,'FontSize',textfontsz)
% title(strcat('Resident = ',num2str(strats(2)),', Invader = ',num2str(invaderstrats(2))),'FontName',fontname,'FontSize',textfontsz);
v=get(ylab,'Position');
set(ylab,'Position',[v(1)+.1 v(2)-.5 v(3)])

ylim=get(gca,'Ylim');
[x,y]=data2norm(0.5,ylim(2));
annotation('textbox',[x-.09 1.05*y .01 .04],'String','(d)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

subplot(3,2,3:4)
nums=cell(Nr,1);
bins=cell(Nr,1);

for i=1:Nr
    [nums{i}, bins{i}]=hist(numrecmat(i,:),0:(max(numrecmat(i,:))+.5));
end

p=zeros(1,Nr);
leglabs=cell(1,Nr);
hold on
for i=1:Nr
    p(i)=bar(bins{i},nums{i}/sum(nums{i}),'FaceColor',mycolormap((i-1)*3+2,:),'EdgeColor','none');
    leglabs{i}=strcat('r = ',num2str(radstohist(i)));
end
xlabel('Number of informed birds','FontName',fontname,'FontSize',textfontsz)
ylab=ylabel('Probability','FontName',fontname,'FontSize',textfontsz);

set(gca,'xlim',[0 N]+.5)
v=get(gca,'Position');
set(gca,'Position',[.1 v(2)+.00 .67 v(4)])
leg=legend(p,leglabs,'FontName',fontname,'FontSize',labfontsz);
legend('boxoff');
v=get(leg,'Position');
set(leg,'Position',[.79 v(2:4)])
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
v=get(ylab,'Position');
set(ylab,'Position',[xlim(1)-.06*diff(xlim) v(2:3)]);

[x,y]=data2norm(.5,ylim(2));
annotation('textbox',[x-.09 1.03*y .01 .04],'String','(b)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','ESSfigure','.pdf');
print(filename,'-dpdf','-r300');
%% figure: prob of invading as fun of # informed and # informed neighbors
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=.5*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);

mycolormap = cbrewer('seq','Reds',50);

subplot(1,2,1)
imagesc(1:maxrec,1:maxrec,transpose(totalprob_invader(1:maxrec,1:maxrec,1)))
v=get(gca,'Position');
set(gca,'Position',[.1 .2 .3 .7])
set(gca,'ydir','normal')
set(gca,'xtick',1:maxrec,'xticklabel',1:maxrec,'ytick',1:maxrec,'yticklabel',1:maxrec,'FontName',fontname,'FontSize',labfontsz)
colormap(mycolormap)
% caxis manual
% caxis([0 1])
ylabel('Number of informed neighbors','FontName',fontname,'FontSize',textfontsz)
xlabel('Number of informed birds in the flock','FontName',fontname,'FontSize',textfontsz)
title(strcat('Resident = ',num2str(strats(1)),', Invader = ',num2str(invaderstrats(1))),'FontName',fontname,'FontSize',textfontsz);

subplot(1,2,2)
imagesc(1:maxrec,1:maxrec,transpose(totalprob_invader(1:maxrec,1:maxrec,2)))
v=get(gca,'Position');
set(gca,'Position',[.5 .2 .3 .7])
set(gca,'ydir','normal')
set(gca,'xtick',1:maxrec,'xticklabel',1:maxrec,'ytick',1:maxrec,'yticklabel',1:maxrec,'FontName',fontname,'FontSize',labfontsz)
colormap(mycolormap)
% caxis manual
% caxis([0 1])
c=colorbar;
v=get(c,'Position');
set(c,'Position',[.85 v(2:4)])
ylim=get(c,'YLabel');
set(ylim,'String','Probability of being eaten','Rotation',270,'Position',[8 .37 1.001]);
ylabel('Number of informed neighbors','FontName',fontname,'FontSize',textfontsz)
xlabel('Number of informed birds in the flock','FontName',fontname,'FontSize',textfontsz)
title(strcat('Resident = ',num2str(strats(2)),', Invader = ',num2str(invaderstrats(2))),'FontName',fontname,'FontSize',textfontsz);

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','mutual_uninvasibility','.pdf');
% print(filename,'-dpdf','-r300');


%% figure: ESS strategies
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=.5*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
fontname='TimesNewRoman';
textfontsz=12;
labfontsz=10;
markersz=3;
lw=2;

cols=cbrewer('seq','Reds',4);

n=size(ESSseaten,2);
k=1;
t=3;

ESSseaten=ESSeaten_total{k};
ESSsgettoeat=ESSgettoeat_total{k};

hold on
for i=1:n
    l=size(ESSseaten{t,i},2);
    p1=plot(radvals(i)*ones(l,1),ESSseaten{t,i},'o','Color',cols(end,:),'MarkerSize',markersz,'LineWidth',lw);
    l=size(ESSsgettoeat{t,i},2);
    if l>0
        p2=plot(radvals(i)*ones(l,1),ESSsgettoeat{t,i},'o','Color','k','MarkerSize',markersz,'LineWidth',lw);
    end
end
xlabel('Radius of signal','FontName',fontname,'FontSize',textfontsz)
ylabel('Number of neighbors','FontName',fontname,'FontSize',textfontsz)

leg=legend([p1 p2],{'Predation','Resources'},'FontName',fontname,'FontSize',labfontsz);
legend('boxoff')
v=get(leg,'Position');
set(leg,'Position',[.15 .4 v(3) v(4)])

set(gca,'FontName',fontname,'FontSize',labfontsz)


set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','ESSneighbors','.pdf');
% print(filename,'-dpdf','-r300');



%% histograms of number of informed birds
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=.5*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
fontname='TimesNewRoman';
textfontsz=12;
labfontsz=10;
markersz=3;
lw=2;

nums=cell(Nr,1);
bins=cell(Nr,1);

cols=cbrewer('seq','Reds',9);

for i=1:Nr
    [nums{i}, bins{i}]=hist(numrecmat(i,:),0:(max(numrecmat(i,:)+1)+.5));
end

p=zeros(1,Nr);
leglabs=cell(1,Nr);
hold on
for i=1:Nr
    p(i)=bar(bins{i},nums{i}/sum(nums{i}),'FaceColor',cols((i-1)*3+1,:));
    leglabs{i}=strcat('r = ',num2str(radstohist(i)));
end
xlabel('Number of informed birds')
ylabel('Probability')
hold off
set(gca,'xlim',[0 N+.5])
leg=legend(p,leglabs);
legend('boxoff');
v=get(leg,'Position');
set(leg,'Position',[.65 v(2:4)])