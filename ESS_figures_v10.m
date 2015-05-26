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
meanlw=3;
thinlw=1;
axislw=.5;

qualcols=cbrewer('qual','Set1',3);
qualcols=qualcols([3 1 2],:);
twocols=cbrewer('qual','Set1',2);
seqcols = cbrewer('seq','Reds',11);
seqcols2 = cbrewer('seq','Reds',7);
seqcols3 = cbrewer('seq','Reds',15);
seqcols4 = cbrewer('seq','Reds',22);
greys=cbrewer('seq','Greys',3);
Ncols=cbrewer('qual','Set1',20);



%% load ESS data

load /Users/eleanorbrush/Desktop/ESS_nummoves=1000_numpermove=1.mat

% load /Users/eleanorbrush/Desktop/greedyopt_simultaneous_T=1_nummoves=1000_numpermove=1_rad=0.2_timesteps=100.mat

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

        receivers = randsample(N,numsigs_permove,'true');
        for j = 1:numsigs_permove

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

resstrats_low = 1:18;
invaderstrats = resstrats_low+1;

[~,costs_to_high] = costs_by_category(N,resstrats_low,invaderstrats,nummoves,radius,b,T);

resstrats_high = 2:19;
invaderstrats = resstrats_high-1;

[costs_to_low,~] = costs_by_category(N,resstrats_high,invaderstrats,nummoves,radius,b,T);
%% find probability of being ranked for low and high strategies
divisions=(2.5:N/5:N);
Nd=length(divisions);
      
low=4;
high=16;
strategy1=low*ones(N,1);
strategy1(1)=high;
scores=signalingevents_returnscores(strategy1,numsigs_permove,nummoves,radius,b,T);
rank=zeros(N,size(scores,2));
for i=1:size(scores,2)
    [~,~,ranked]=unique(scores(:,i));
    rank(:,i)=ranked;
end
[lowres_low_count,~]=hist(rank(2,:),divisions);
[lowres_high_count,~]=hist(rank(1,:),divisions);

strategy2=high*ones(N,1);
strategy2(1)=low;
scores=signalingevents_returnscores(strategy2,numsigs_permove,nummoves,radius,b,T);
rank=zeros(N,size(scores,2));
for i=1:size(scores,2)
    [~,~,ranked]=unique(scores(:,i));
    rank(:,i)=ranked;
end
[highres_low_count,~]=hist(rank(1,:),divisions);
[highres_high_count,~]=hist(rank(2,:),divisions);

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
t=1;

subplot(2,2,1)
hold on

set(gca,'xlim',[-.05 1.45]);
xlim=get(gca,'xlim');
% xdiv=[.25 1.05];
% plot(xdiv(1)*ones(2,1),ylim,'--k','LineWidth',1)
% plot(xdiv(2)*ones(2,1),ylim,'--k','LineWidth',1)
% fill([xlim(1) xdiv(1) xdiv(1) xlim(1)],[0 0 N N],qualcols(1,:),'EdgeColor','none')
% fill([xdiv(1) xdiv(2) xdiv(2) xdiv(1)],[0 0 N N],qualcols(2,:),'EdgeColor','none')
% fill([xdiv(2) xlim(2) xlim(2) xdiv(2)],[0 0 N N],qualcols(3,:),'EdgeColor','none')
for i=1:n
    l=size(ESSseaten{t,i},2);
    if l>0
        p1=plot(radvals(i)*ones(l,1),ESSseaten{t,i},'o','Color','k','MarkerFaceColor','k','MarkerSize',markersz,'LineWidth',lw);
    end
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
set(gca,'xlim',[0 1])
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
colormap(flipud(seqcols))
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
set(gca,'FontName',fontname,'FontSize',labfontsz)

apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String','(b)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

subplot(2,2,3)
hold on
barcols=cbrewer('div','RdBu',4);

width=.2;
hold on
bar((1:Nd)-3*width/2,lowres_low_count/nummoves,width,'FaceColor',barcols(1,:),'EdgeColor','k')
bar((1:Nd)-width/2,lowres_high_count/nummoves,width,'FaceColor',barcols(4,:),'EdgeColor','k')
bar((1:Nd)+width/2,highres_low_count/nummoves,width,'FaceColor',barcols(2,:),'EdgeColor','k')
bar((1:Nd)+3*width/2,highres_high_count/nummoves,width,'FaceColor',barcols(3,:),'EdgeColor','k')
set(gca,'ylim',[0 1])
[leg,hobj]=legend({'Low','High','Low resident','High resident'},'FontName',fontname,'FontSize',labfontsz);
legend 'boxoff'
v=get(leg,'position');
set(leg,'position',[v(1)-.112 v(2)+.07 v(3) v(4)])
legend('boxoff')
textobj=findobj(hobj,'type','patch');
graycols=cbrewer('div','RdGy',4);
graycols=graycols([1 2 4 3],:);
for i=1:2
    set(textobj(i),'xdata',[.29 .29 .36 .36 .29])
end
for i=3:4
    set(textobj(i),'xdata',[.29 .29 .36 .36 .29],'FaceColor',graycols(i,:))
end
xlabel('Rank','FontName',fontname,'FontSize',textfontsz)
ylabel('Probability','FontName',fontname,'FontSize',textfontsz)
set(gca,'FontName',fontname,'FontSize',labfontsz)
set(gca,'xtick',[1:5],'XTickLabel',{'1-4','5-8','9-12','13-16','17-20'})

apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String','(c)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

subplot(2,2,4)
hold on
mybar=bar(costs_to_high([1 4 5],:)','stack');
for k=1:3
  set(mybar(k),'facecolor',qualcols(k,:))
end
xlim=get(gca,'xlim');
plot(xlim,1/N*ones(2,1),'--k','LineWidth',axislw)
[leg,hobj]=legend(mybar, {'Neither','Both','Multiple'},'FontName',fontname,'FontSize',labfontsz);
v=get(leg,'position');
set(leg,'position',[v(1)-.2 v(2)+.05 v(3) v(4)])
legend('boxoff')
textobj=findobj(hobj,'type','patch');
for i=1:3
    set(textobj(i),'xdata',[.3 .3 .4 .4 .3])
end
set(gca,'ylim',[0 1])
xlabel('Resident strategy','FontName',fontname,'FontSize',textfontsz)
ylabel('Probability','FontName',fontname,'FontSize',textfontsz)
set(gca,'xtick',1:2:length(resstrats_low),'xticklabel',min(resstrats_low):2:max(resstrats_low))
set(gca,'FontName',fontname,'FontSize',labfontsz)

apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String','(d)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')




set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','ESSfigure_full','.pdf');
print(filename,'-dpdf','-r300');

%% figure: ESS strategies as a function of radius
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=.25*w;
set(gcf,'Units','inches');
set(gcf,'Position',[11.5 3 w h]);

xoffset=.1;
yoffset=.06;


n=size(ESSseaten,2);
k=1;
t=1;


subplot(1,2,1)
r=2;
fitnesseaten=storefitnesseaten{t,r};
imagesc(strats,strats,transpose(fitnesseaten))
set(gca,'xtick',[1 5:5:15 N-1],'ytick',[1 5:5:15 N-1])
set(gca,'ydir','normal')

myseqmap=cbrewer('seq', 'YlOrRd',9);
myseqmap=myseqmap(end:-1:1,:);
colormap(flipud(seqcols))
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
set(gca,'FontName',fontname,'FontSize',labfontsz)

apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String','(a)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

subplot(1,2,2)
hold on

set(gca,'xlim',[-.05 1.45]);
xlim=get(gca,'xlim');
for i=1:n
    l=size(ESSseaten{t,i},2);
    if l>0
        p1=plot(radvals(i)*ones(l,1),ESSseaten{t,i},'o','Color','k','MarkerFaceColor','k','MarkerSize',markersz,'LineWidth',lw);
    end
    l=size(ESSsgettoeat{t,i},2);
    if l>0
        p2=plot(radvals(i)*ones(l,1),ESSsgettoeat{t,i},'s','Color','k','MarkerFaceColor','k','MarkerSize',markersz,'LineWidth',lw);
    end
    l=size(ESSsboth{t,i},2);
    if l>0
        p3=plot(radvals(i)*ones(l,1),ESSsboth{t,i},'d','Color','k','MarkerFaceColor','k','MarkerSize',markersz,'LineWidth',lw);
    end
end
xlab=xlabel('Radius of signal','FontName',fontname,'FontSize',textfontsz);
ylab=ylabel('ESS','FontName',fontname,'FontSize',textfontsz);
v=get(ylab,'Position');
set(ylab,'Position',[xlim(1)-.05*diff(xlim) v(2:3)]);
set(gca,'xlim',[0 1])
ylim=get(gca,'Ylim');
plot(xlim(1)*ones(2,1),ylim,'k','LineWidth',axislw);
plot(xlim,(ylim(1)-.0)*ones(2,1),'k','LineWidth',axislw);

leg=legend([p1 p2 p3],{'Predation','Resources','Both'},'FontName',fontname,'FontSize',labfontsz);
legend('boxoff')
legpos=get(leg,'Position');
set(leg,'Position',[.82 .47 legpos(3) legpos(4)])

set(gca,'FontName',fontname,'FontSize',labfontsz,'tickdir','out')
v=get(gca,'Position');
% set(gca,'Position',[.6 v(2:3) .8])
% set(gca,'Position',[.6 .25 .3347 .6242])
apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset+.02 apos(2)+apos(4)+yoffset .01 .04],'String','(b)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','ESSfigure','.pdf');
% print(filename,'-dpdf','-r300');

%% how ESS depends on robustness
figure
set(gcf,'Color','w')
w=6.83;
h=1/3*w;
set(gcf,'Units','inches');
set(gcf,'Position',[11.5 3 w h]);
subplot(1,2,1)

imagesc(radvals,homogenstrats,rhos_forced)
caxis manual
set(gca,'ydir','normal')

colormap(seqcols)
cb=colorbar;

ylab=ylabel(cb,'Robustness','FontName',fontname,'FontSize',labfontsz);
set(ylab,'Rotation',270);
v=get(ylab,'Position');
set(ylab,'Position',[v(1)+.1 v(2) v(3)])
firstxlab=xlabel('Radius of signal','FontName',fontname,'FontSize',textfontsz);
ylabel('Strategy','FontName',fontname,'FontSize',textfontsz);

set(gca,'FontName',fontname,'FontSize',labfontsz)

apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String','(a)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

v=get(gca,'Position');
set(gca,'Position',[.08 .18 .27 .72])


subplot(1,2,2)
hold on
t=1;
n=size(ESSseaten,2);
for i=1:12
    l=length(ESSseaten{t,i});
    plot(corrlengths_forced(12,i)*ones(l,1),ESSseaten{t,i},'ok','MarkerFaceColor','k','LineWidth',lw,'MarkerSize',markersz)
%     plot(corrlengths_forced(ESSseaten{t,i},i),ESSseaten{t,i},'ok','MarkerFaceColor','k','LineWidth',lw,'MarkerSize',markersz)
end
set(gca,'ylim',[13 20])
xlabel('Robustness','FontName',fontname,'FontSize',textfontsz)
ylabel('ESS','FontName',fontname,'FontSize',textfontsz)

set(gca,'FontName',fontname,'FontSize',labfontsz)

apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String','(b)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')


set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename='/Users/eleanorbrush/Desktop/robustness_vs_ESS.pdf';
% print(filename,'-dpdf','-r300')

%% H2 / correlation length relationship
figure
set(gcf,'Color','w')
w=6.83;
h=1/2*w;
set(gcf,'Units','inches');
set(gcf,'Position',[11.5 3 w h]);

xoffset=.12;
yoffset=.06;

stratvals=1:19;

subplot(2,2,1)
imagesc(radvals,stratvals,rhos_forced)
set(gca,'ydir','normal')
set(gca,'FontName',fontname,'FontSize',labfontsz)
caxis manual
colormap(seqcols)
cb=colorbar;
set(cb,'FontName',fontname,'FontSize',labfontsz)
ylab=ylabel(cb,'Robustness','FontName',fontname,'FontSize',labfontsz);
set(ylab,'Rotation',270);
v=get(ylab,'Position');
set(ylab,'Position',[v(1)+.1 v(2) v(3)])

xlabel('Radius','FontName',fontname,'FontSize',textfontsz)
ylabel('Strategy','FontName',fontname,'FontSize',textfontsz)

apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String','(a)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

subplot(2,2,2)
imagesc(radvals,stratvals,corrlengths_forced)
set(gca,'ydir','normal')
set(gca,'FontName',fontname,'FontSize',labfontsz)
caxis manual
colormap(seqcols)
cb=colorbar;
set(cb,'FontName',fontname,'FontSize',labfontsz)
ylab=ylabel(cb,'Correlation length','FontName',fontname,'FontSize',labfontsz);
set(ylab,'Rotation',270);
v=get(ylab,'Position');
set(ylab,'Position',[v(1)+.1 v(2) v(3)])

xlabel('Radius','FontName',fontname,'FontSize',textfontsz)
ylabel('Strategy','FontName',fontname,'FontSize',textfontsz)

apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String','(b)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')


subplot(2,2,3:4)
hold on

% toplot=1:2:length(0:.1:1);
toplot=1:2:14;

k=6;
for i=1:length(toplot)
    if i<8
    k=6;
    else k=1;
    end
    plot(1./H2s_forced(k:end,toplot(i)),corrlengths_forced(k:end,toplot(i)),'-o','Color',seqcols4(20-i,:),'MarkerSize',markersz,'LineWidth',lw)
end

leglabs=cell(length(toplot),1);
for i=1:length(toplot)
    leglabs{i}=num2str(radvals(toplot(i)));
end
set(gca,'FontName',fontname,'FontSize',labfontsz)
% 
% toplot=2:18;
% 
% for i=1:length(toplot)
%     plot(rhos_forced(toplot(i),:),corrlengths_forced(toplot(i),:)/max(corrlengths_forced(toplot(i),:)),'-o','Color',seqcols4(21-i,:),'MarkerSize',markersz,'LineWidth',lw)
% end
% 
% 
% % set(gca,'ylim',[0 .5])
% set(gca,'FontName',fontname,'FontSize',labfontsz)
% 
% box off
% 
% leglabs=cell(length(toplot),1);
% for i=1:length(toplot)
%     leglabs{i}=num2str(stratvals(toplot(i)));
% end

[leg,hobj]=legend(leglabs);
legend('boxoff')
legpos=get(leg,'position');
set(leg,'position',[.83 .17 legpos(3:4)],'FontName',fontname,'FontSize',labfontsz)

textobj = findobj(hobj, 'type', 'line');
for i=1:2:length(textobj)
    set(textobj(i),'XData',[.2 .4])
end

apos=get(gca,'Position');
set(gca,'Position',[apos(1) .15 apos(3:4)])
apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String','(c)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

xlabel('Robustness','FontName',fontname,'FontSize',textfontsz)
ylabel('Correlation length','FontName',fontname,'FontSize',textfontsz)

axpos=get(gca,'position');
% set(gca,'position',[.08 axpos(2:4)])


set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename='/Users/eleanorbrush/Desktop/H2_v_corrlength.pdf';
print(filename,'-dpdf','-r300')

%% H2 / correlation length heat maps
figure
set(gcf,'Color','w')
w=6.83;
h=1/4*w;
set(gcf,'Units','inches');
set(gcf,'Position',[11.5 3 w h]);

xoffset=.08;
yoffset=.06;

stratvals=1:19;

subplot(1,2,1)
imagesc(radvals,stratvals,rhos_forced)
set(gca,'ydir','normal')
set(gca,'FontName',fontname,'FontSize',labfontsz)
set(gca,'ytick',[1 10 19],'yticklabel',[1 10 19])
caxis manual
colormap(seqcols)
cb=colorbar;
set(cb,'FontName',fontname,'FontSize',labfontsz)
ylab=ylabel(cb,'Robustness','FontName',fontname,'FontSize',labfontsz);
set(ylab,'Rotation',270);
v=get(ylab,'Position');
set(ylab,'Position',[v(1)+.1 v(2) v(3)])

xlabel('Radius','FontName',fontname,'FontSize',textfontsz)
ylabel('Strategy','FontName',fontname,'FontSize',textfontsz)

apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String','(a)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

subplot(1,2,2)
imagesc(radvals,stratvals,corrlengths_forced)
set(gca,'ydir','normal')
set(gca,'FontName',fontname,'FontSize',labfontsz)
set(gca,'ytick',[1 10 19],'yticklabel',[1 10 19])
caxis manual
colormap(seqcols)
cb=colorbar;
set(cb,'FontName',fontname,'FontSize',labfontsz)
ylab=ylabel(cb,'Correlation length','FontName',fontname,'FontSize',labfontsz);
set(ylab,'Rotation',270);
v=get(ylab,'Position');
set(ylab,'Position',[v(1)+.1 v(2) v(3)])

xlabel('Radius','FontName',fontname,'FontSize',textfontsz)
ylabel('Strategy','FontName',fontname,'FontSize',textfontsz)

apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String','(b)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')


set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename='/Users/eleanorbrush/Desktop/H2_and_corrlength.pdf';
% print(filename,'-dpdf','-r300')

%% fitness heatmaps
j=j+1;i=1;

subplot(1,4,1)
fitness=storefitnesseaten{i,j};
imagesc(strats,strats,transpose(log(fitness)));set(gca,'ydir','normal');
m=min(min(log(fitness)));M=max(max(log(fitness)));
M=max(abs(m),M);caxis manual;caxis([-M M])
colorbar

subplot(1,4,2)
fitness=storefitnessgettoeat{i,j};
imagesc(strats,strats,transpose(log(fitness)));set(gca,'ydir','normal');
m=min(min(log(fitness)));M=max(max(log(fitness)));
M=max(abs(m),M);caxis manual;caxis([-M M])
colorbar

subplot(1,4,3)
fitness=storefitnessboth{i,j};
imagesc(strats,strats,transpose(log(fitness)));set(gca,'ydir','normal');
m=min(min(log(fitness)));M=max(max(log(fitness)));
M=max(abs(m),M);caxis manual;caxis([-M M])
colorbar

subplot(1,4,4)
fitness=storefitnessgenerous{i,j};
imagesc(strats,strats,transpose(log(fitness)));set(gca,'ydir','normal');
m=min(min(log(fitness)));M=max(max(log(fitness)));
M=max(abs(m),M);caxis manual;caxis([-M M])
colorbar
