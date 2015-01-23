%% setup
N = 20;
nummoves = 1000;
numsigs_permove = 1;
b = 1;

fontname='Times';
textfontsz=12;
labfontsz=9;

markersz=3;
lw=2;
axislw=.5;
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
%% find number of receivers for different radius
radstohist = [.1 .5 1.2];
T = 1; %#ok<NASGU>
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
%% find costs to low and high strategies
N=20;
radius = .1;
nummoves=1000;
T = 1;
% strats=18*ones(17,1);
mainstrats = 3:17;
% strats=[14:17];
% invaderstrats = 1:17;
% invaderstrats = 18*ones(Ls,1);
invaderstrats = mainstrats+1;
% invaderstrats = fliplr(strats);

[~,costs_to_high] = costs_by_category(N,mainstrats,invaderstrats,nummoves,radius,b,T);

% strats=18*ones(17,1);
mainstrats = 3:17;
% strats=[14:17];
% invaderstrats = 1:17;
% invaderstrats = 18*ones(Ls,1);
invaderstrats = mainstrats-1;
% invaderstrats = fliplr(strats);

[costs_to_low,~] = costs_by_category(N,mainstrats,invaderstrats,nummoves,radius,b,T);


%% figure: ESS strategies and prob eaten comparison together
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=.5*w;
set(gcf,'Units','inches');
set(gcf,'Position',[11.5 3 w h]);

xoffset=.1;
yoffset=.06;

n=size(ESSseaten,2);
k=1;
t=3;

mycolormap = cbrewer('seq','Reds',11);

ESSseaten=ESSeaten_total{k};
ESSsgettoeat=ESSgettoeat_total{k};

subplot(2,2,1)
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
    p1=plot(radvals(i)*ones(l,1),ESSseaten{t,i},'o','Color','k','MarkerFaceColor','k','MarkerSize',markersz,'LineWidth',lw);
    l=size(ESSsgettoeat{t,i},2);
    if l>0
        p2=plot(radvals(i)*ones(l,1),ESSsgettoeat{t,i},'s','Color','k','MarkerFaceColor','k','MarkerSize',markersz,'LineWidth',lw);
    end
end
xlab=xlabel('Radius of signal','FontName',fontname,'FontSize',textfontsz);
ylab=ylabel('ESS','FontName',fontname,'FontSize',textfontsz);
v=get(xlab,'Position');
set(xlab,'Position',[v(1) -4.5 v(3)]);
v=get(ylab,'Position');
set(ylab,'Position',[xlim(1)-.1*diff(xlim) v(2:3)]);
ylim=get(gca,'Ylim');
plot(xlim(1)*ones(2,1),ylim,'k','LineWidth',axislw);
plot(xlim,(ylim(1)-.0)*ones(2,1),'k','LineWidth',axislw);

leg=legend([p1 p2],{'Predation','Resources'},'FontName',fontname,'FontSize',labfontsz);
legend('boxoff')
legpos=get(leg,'Position');
set(leg,'Position',[.18 .73 legpos(3) legpos(4)])

set(gca,'FontName',fontname,'FontSize',labfontsz,'tickdir','out')

apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String','(a)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

subplot(2,2,2)
r=2;
fitnesseaten=storefitnesseaten{t,r};
imagesc(strats,strats,transpose(fitnesseaten))
set(gca,'xtick',[1 5:5:15 N-1],'ytick',[1 5:5:15 N-1])
set(gca,'ydir','normal')

myseqmap=cbrewer('seq', 'YlOrRd',9);
myseqmap=myseqmap(end:-1:1,:);
colormap(flipud(mycolormap))
cb=colorbar;
caxis manual
caxis([0 1])
yl=ylabel(cb,'Relative fitness','FontName',fontname,'FontSize',textfontsz);
v=get(yl,'Position');
set(yl,'Rotation',270,'Position',[v(1)+.5 v(2:3)])
hold on
plot(ESSseaten{t,r},ESSseaten{t,r},'o','Color','k','MarkerFaceColor','k','MarkerSize',markersz,'LineWidth',lw);
xlabel('Resident strategy','FontName',fontname,'FontSize',textfontsz)
ylabel('Invader strategy','FontName',fontname,'FontSize',textfontsz)
v=get(gca,'Position');
set(gca,'FontName',fontname,'FontSize',labfontsz,'Position',[v(1)+.05 v(2:4)])

apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String','(b)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')


subplot(2,2,3)
hold on
plot(mainstrats,costs_to_low(3,:),'-o','Color',mycolormap(end-3,:),'LineWidth',lw)
plot(mainstrats,costs_to_high(end,:),'-ok','LineWidth',lw)
set(gca,'xlim',[mainstrats(1)-1 mainstrats(end)+1])
xlabel('Resident strategy','FontName',fontname,'FontSize',textfontsz)
ylabel('Predation risk','FontName',fontname,'FontSize',textfontsz)
set(gca,'FontName',fontname,'FontSize',labfontsz)
[leg,hobj]=legend({'Lower','Higher'},'FontName',fontname,'FontSize',labfontsz);
legend('boxoff')

ylim=get(gca,'ylim');
[x,y]=data2norm(3,mean(ylim));
v=get(leg,'position');
set(leg,'position',[x y v(3) v(4)]);
textobj = findobj(hobj, 'type', 'line');
for i=[1 3]
    set(textobj(i),'XData',[.2 .45])
end

apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String','(c)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

subplot(2,2,4)
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
    leglabs{i}=strcat('r =  ',num2str(radstohist(i)));
end
xlabel('Number of informed birds','FontName',fontname,'FontSize',textfontsz)
ylab=ylabel('Probability','FontName',fontname,'FontSize',textfontsz);
set(gca,'xlim',[0 N]+.5)

[leg,hobj]=legend(p,leglabs,'FontName',fontname,'FontSize',labfontsz);
legend('boxoff');
v=get(leg,'Position');
set(leg,'Position',[.74 .23 v(3:4)])
textobj=findobj(hobj,'type','patch');
for i=1:3
    set(textobj(i),'xdata',[.3 .3 .4801 .4801 .3])
end

xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
plot(xlim(1)*ones(2,1),ylim,'k','LineWidth',axislw);
plot(xlim,(ylim(1)-.0)*ones(2,1),'k','LineWidth',axislw);

set(gca,'FontName',fontname,'FontSize',labfontsz,'tickdir','out')
v=get(gca,'Position');
set(gca,'FontName',fontname,'FontSize',labfontsz,'tickdir','out','Position',[v(1)+.05 v(2:4)])
apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String','(d)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')


set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','ESSfigure','.pdf');
print(filename,'-dpdf','-r300');
