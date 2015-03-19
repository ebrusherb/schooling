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

qualcols=cbrewer('qual','Set1',5);
qualcols=qualcols([1 3 2 4 5],:);
seqcols = cbrewer('seq','Reds',11);
seqcols2 = cbrewer('seq','Reds',7);


xoffset=.1;
yoffset=.06;

%% load ESS data

load /Users/eleanorbrush/Desktop/ESS_nummoves=1000_numpermove=1.mat

load /Users/eleanorbrush/Desktop/greedyopt_simultaneous_T=1_nummoves=1000_numpermove=1_rad=0.2_timesteps=100.mat

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


%% figure: ESS strategies and prob eaten comparison together
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=.5*w;
set(gcf,'Units','inches');
set(gcf,'Position',[11.5 3 w h]);

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


subplot(2,2,4)
hold on
mybar=bar(costs_to_high','stack');
for k=1:5
  set(mybar(k),'facecolor',qualcols(k,:))
end
xlim=get(gca,'xlim');
plot(xlim,1/N*ones(2,1),'--k','LineWidth',axislw)
[leg,hobj]=legend(mybar, {'Neither','Resident','Invader','Both','Multiple'},'FontName',fontname,'FontSize',labfontsz);
v=get(leg,'position');
set(leg,'position',[v(1)-.2 v(2)+.05 v(3) v(4)])
legend('boxoff')
textobj=findobj(hobj,'type','patch');
for i=1:5
    set(textobj(i),'xdata',[.3 .3 .4 .4 .3])
end
set(gca,'ylim',[0 1])
xlabel('Resident strategy','FontName',fontname,'FontSize',textfontsz)
ylabel('Probability','FontName',fontname,'FontSize',textfontsz)
set(gca,'xtick',1:2:length(resstrats_low),'xticklabel',min(resstrats_low):2:max(resstrats_low))
set(gca,'FontName',fontname,'FontSize',labfontsz)

apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String','(d)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

subplot(2,2,3)
hold on
mybar=bar(costs_to_low([1 3 2 4 5],:)','stack');
for k=1:5
  set(mybar(k),'facecolor',qualcols(k,:))
end
[leg,hobj]=legend(mybar, {'Neither','Resident','Invader','Both','Multiple'},'FontName',fontname,'FontSize',labfontsz);
v=get(leg,'position');
set(leg,'position',[v(1)-.2 v(2)+.05 v(3) v(4)])
legend('boxoff')
textobj=findobj(hobj,'type','patch');
for i=1:5
    set(textobj(i),'xdata',[.3 .3 .4 .4 .3])
end
xlim=get(gca,'xlim');
plot(xlim,1/N*ones(2,1),'--k','LineWidth',axislw)
set(gca,'ylim',[0 1])
xlabel('Resident strategy','FontName',fontname,'FontSize',textfontsz)
ylabel('Probability','FontName',fontname,'FontSize',textfontsz)
set(gca,'xtick',1:2:length(resstrats_high),'xticklabel',min(resstrats_high):2:max(resstrats_high))
set(gca,'FontName',fontname,'FontSize',labfontsz)

apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String','(c)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')


set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','ESSfigure','.pdf');
print(filename,'-dpdf','-r300');


%% how ESS depends on robustness
figure
set(gcf,'Color','w')
w=6.83;
h=1/3*w;
set(gcf,'Units','inches');
set(gcf,'Position',[11.5 3 w h]);
subplot(1,2,1)
% toplot=corrlengths_forced;
% toplot=corrlengths./repmat([1 radvals(2:end)],Lh,1);
toplot=ones(14,15)./(groupconsensus_forced.*power(repmat(col(6:19),1,15),1/2));
toplot_vec=1./(groupconsensus_forced(1,:)*sqrt(6));

imagesc(radvals,homogenstrats,toplot)
caxis manual
set(gca,'ydir','normal')

% a=min(min(toplot))-1;
A=max(max(toplot))-1;
% A=max(abs(a),abs(A))+1;
% caxis([0 max(A+1,.5)])
% caxis([-A A])

colormap(seqcols)
cb=colorbar;

ylab=ylabel(cb,'Normalized robustness','FontName',fontname,'FontSize',labfontsz);
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
    plot(toplot_vec(i)*ones(l,1),ESSseaten{t,i},'ok','MarkerFaceColor','k','LineWidth',lw,'MarkerSize',markersz)
end
set(gca,'ylim',[13 20])
xlabel('Normalized robustness','FontName',fontname,'FontSize',textfontsz)
ylabel('ESS','FontName',fontname,'FontSize',textfontsz)

set(gca,'FontName',fontname,'FontSize',labfontsz)

apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String','(b)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')


set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename='/Users/eleanorbrush/Desktop/robustness_vs_ESS.pdf';
print(filename,'-dpdf','-r300')

%% H2 / correlation length relationship
figure
set(gcf,'Color','w')
w=6.83;
h=1/2*w;
set(gcf,'Units','inches');
set(gcf,'Position',[11.5 3 w h]);

subplot(2,2,1)
imagesc(radvals,6:19,ones(14,15)./(groupconsensus_forced.*power(repmat(col(6:19),1,15),1/2)))
set(gca,'ydir','normal')
set(gca,'FontName',fontname,'FontSize',labfontsz)
caxis manual
colormap(seqcols)
cb=colorbar;
set(cb,'FontName',fontname,'FontSize',labfontsz)
ylab=ylabel(cb,'Normalized robustness','FontName',fontname,'FontSize',labfontsz);
set(ylab,'Rotation',270);
v=get(ylab,'Position');
set(ylab,'Position',[v(1)+.1 v(2) v(3)])

xlabel('Radius','FontName',fontname,'FontSize',textfontsz)
ylabel('Strategy','FontName',fontname,'FontSize',textfontsz)

apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String','(a)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

subplot(2,2,2)
imagesc(radvals,6:19,corrlengths_forced)
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
radvals=0:.1:1.4;
toplot=1:2:length(0:.1:1);

% plot(ones(14,1)./(groupconsensus(:,1).*col(6:19)),corrlengths(:,1),'-o','Color',seqcols2(end,:),'MarkerSize',markersz,'LineWidth',lw)

for i=1:length(toplot)
    plot(ones(14,1)./(groupconsensus_forced(:,toplot(i)).*sqrt(col(6:19))),corrlengths_forced(:,toplot(i)),'-o','Color',seqcols2(8-i,:),'MarkerSize',markersz,'LineWidth',lw)
end
set(gca,'ylim',[0 .5])
set(gca,'FontName',fontname,'FontSize',labfontsz)

box off

leglabs=cell(length(toplot),1);
% leglabs{1}='No signal';
for i=1:length(toplot)
    leglabs{i}=num2str(radvals(toplot(i)));
end

[leg,hobj]=legend(leglabs);
legend('boxoff')
legpos=get(leg,'position');
set(leg,'position',[.83 .17 legpos(3:4)],'FontName',fontname,'FontSize',labfontsz)

textobj = findobj(hobj, 'type', 'line');
for i=[1:2:length(textobj)]
    set(textobj(i),'XData',[.2 .4])
end

apos=get(gca,'Position');
set(gca,'Position',[apos(1) .15 apos(3:4)])
apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String','(c)','FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

xlabel('Normalized robustness','FontName',fontname,'FontSize',textfontsz)
ylabel('Correlation length','FontName',fontname,'FontSize',textfontsz)

axpos=get(gca,'position');
% set(gca,'position',[.08 axpos(2:4)])


set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename='/Users/eleanorbrush/Desktop/H2_v_corrlength.pdf';
print(filename,'-dpdf','-r300')