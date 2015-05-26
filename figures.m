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
Ncols=cbrewer('qual','Set1',N);
%% figure: H2 / correlation length heat maps
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
%% figure: greedy opt and group props over time
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=4/3*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
fontname='TimesNewRoman';
textfontsz=12;
markersz=3;

xoffset=.08;
yoffset=-.02;

letters={'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)','(j)','(k)','(l)'};

timesteps=size(evolution_eaten,3);
M=max([max(maxt_eaten),max(maxt_gettoeat),max(maxt_both),max(maxt_generous)]);
M=min(timesteps,M+5);

rhomax=max(max(max(rhos_eaten)),max(max(rhos_gettoeat)))+.03;
corrlengthmax=max(max(max(corrlengths_eaten)),max(max(corrlengths_eaten)))+.08;

subplot(4,3,1)
i=2;
hold on
for k=1:N
    plot(1:timesteps,col(evolution_eaten(i,k,:)),'LineWidth',thinlw,'Color',Ncols(k,:))
end
set(gca,'xlim',[0 M])
box off
set(gca,'xtick',0:25:M,'ytick',0:5:N,'tickdir','out')
set(gca,'FontName',fontname,'FontSize',textfontsz)
xlabel('Timesteps','FontName',fontname,'FontSize',textfontsz)
ylabel('Number of neighbors','FontName',fontname,'FontSize',textfontsz)
apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String',letters{1},'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

        
subplot(4,3,2)
plot(1:timesteps,rhos_eaten,'Color',greys(2,:),'LineWidth',thinlw)
box off
hold on

plot(1:timesteps,mean(rhos_eaten),'Color','k','LineWidth',lw)
set(gca,'xlim',[0 M],'ylim',[0 rhomax+.005])
set(gca,'xtick',0:25:M,'ytick',0:.05:.15,'tickdir','out')
set(gca,'FontName',fontname,'FontSize',textfontsz)
xlabel('Timesteps','FontName',fontname,'FontSize',textfontsz)
ylabel('Robustness','FontName',fontname,'FontSize',textfontsz)
apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String',letters{2},'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

subplot(4,3,3)
plot(1:timesteps,corrlengths_eaten,'Color',greys(2,:),'LineWidth',thinlw)
box off
hold on
plot(1:timesteps,mean(corrlengths_eaten),'Color','k','LineWidth',lw)
set(gca,'xlim',[0 M],'ylim',[0 corrlengthmax])
set(gca,'xtick',0:25:M,'ytick',0:.2:.4,'tickdir','out')
set(gca,'FontName',fontname,'FontSize',textfontsz)
xlabel('Timesteps','FontName',fontname,'FontSize',textfontsz)
ylabel('Correlation length','FontName',fontname,'FontSize',textfontsz)
apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String',letters{3},'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

subplot(4,3,4)
hold on
for k=1:N
    plot(1:timesteps,col(evolution_gettoeat(1,k,:)),'LineWidth',thinlw,'Color',Ncols(k,:))
end
box off
set(gca,'xtick',0:25:M,'ytick',0:5:N,'tickdir','out')
set(gca,'xlim',[0 M])
set(gca,'FontName',fontname,'FontSize',textfontsz)
apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String',letters{4},'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

 
subplot(4,3,5)
plot(1:timesteps,rhos_gettoeat,'Color',greys(2,:),'LineWidth',thinlw)
hold on
plot(1:timesteps,mean(rhos_gettoeat),'Color','k','LineWidth',lw)
box off
set(gca,'xlim',[0 M],'ylim',[0 rhomax+.005])
set(gca,'xtick',0:25:M,'ytick',0:.05:.15,'tickdir','out')
set(gca,'FontName',fontname,'FontSize',textfontsz)
apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String',letters{5},'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')
% annotation('textbox',[apos(1)+.5*apos(3) apos(2)+apos(4)+yoffset .01 .04],'String','Resources','FitBoxtoText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

subplot(4,3,6)
plot(1:timesteps,corrlengths_gettoeat,'Color',greys(2,:),'LineWidth',thinlw)
box off
hold on
plot(1:timesteps,mean(corrlengths_gettoeat),'Color','k','LineWidth',lw)
box off
set(gca,'xlim',[0 M],'ylim',[0 corrlengthmax])
set(gca,'xtick',0:25:M,'ytick',0:.2:.4,'tickdir','out')
set(gca,'FontName',fontname,'FontSize',textfontsz)
apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String',letters{6},'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

i=9;
subplot(4,3,7)
hold on
for k=1:N
    plot(1:timesteps,col(evolution_both(i,k,:)),'LineWidth',thinlw,'Color',Ncols(k,:))
end
box off
set(gca,'xlim',[0 M])
set(gca,'xtick',0:25:M,'ytick',0:5:N,'tickdir','out')
apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String',letters{7},'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')


set(gca,'FontName',fontname,'FontSize',textfontsz)
subplot(4,3,8)
plot(1:timesteps,rhos_both,'Color',greys(2,:),'LineWidth',thinlw)
hold on
plot(1:timesteps,mean(rhos_both),'Color','k','LineWidth',lw)
box off
set(gca,'xlim',[0 M],'ylim',[0 rhomax+.005])
set(gca,'xtick',0:25:M,'ytick',0:.05:.15,'tickdir','out')
set(gca,'FontName',fontname,'FontSize',textfontsz)
apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String',letters{8},'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

subplot(4,3,9)
plot(1:timesteps,corrlengths_both,'Color',greys(2,:),'LineWidth',thinlw)
box off
hold on
plot(1:timesteps,mean(corrlengths_both),'Color','k','LineWidth',lw)
box off
set(gca,'xlim',[0 M],'ylim',[0 corrlengthmax])
set(gca,'xtick',0:25:M,'ytick',0:.2:.4,'tickdir','out')
set(gca,'FontName',fontname,'FontSize',textfontsz)
apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String',letters{9},'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

i=10;
subplot(4,3,10)
hold on
for k=1:N
    plot(1:timesteps,col(evolution_generous(i,k,:)),'LineWidth',thinlw,'Color',Ncols(k,:))
end
box off
set(gca,'xlim',[0 M])
set(gca,'xtick',0:25:M,'ytick',0:5:N,'tickdir','out')
apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String',letters{10},'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')


set(gca,'FontName',fontname,'FontSize',textfontsz)
subplot(4,3,11)
plot(1:timesteps,rhos_generous,'Color',greys(2,:),'LineWidth',thinlw)
hold on
plot(1:timesteps,mean(rhos_generous),'Color','k','LineWidth',lw)
box off
set(gca,'xlim',[0 M],'ylim',[0 rhomax+.005])
set(gca,'xtick',0:25:M,'ytick',0:.05:.15,'tickdir','out')
set(gca,'FontName',fontname,'FontSize',textfontsz)
apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String',letters{11},'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')

subplot(4,3,12)
plot(1:timesteps,corrlengths_generous,'Color',greys(2,:),'LineWidth',thinlw)
box off
hold on
plot(1:timesteps,mean(corrlengths_generous),'Color','k','LineWidth',lw)
box off
set(gca,'xlim',[0 M],'ylim',[0 corrlengthmax])
set(gca,'xtick',0:25:M,'ytick',0:.2:.4,'tickdir','out')
set(gca,'FontName',fontname,'FontSize',textfontsz)
apos=get(gca,'Position');
annotation('textbox',[apos(1)-xoffset apos(2)+apos(4)+yoffset .01 .04],'String',letters{12},'FitBoxToText','on','FontSize',textfontsz,'FontName',fontname,'EdgeColor','none','VerticalAlignment','middle','HorizontalAlignment','left')


set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','greedyoptneighbors_radius=',num2str(radius*10),'.pdf');
% print(filename,'-dpdf','-r300');

print('/Users/eleanorbrush/Desktop/test.pdf','-dpdf','-r300')
%% figure: hypergeometric ESS
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


n=size(ESSseaten_hyper,2);
k=1;
t=1;


subplot(1,2,1)
r=2;
fitnesseaten_hyper=storefitnesseaten_hyper{t,r};
imagesc(strats,strats,transpose(fitnesseaten_hyper))
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
plot(ESSseaten_hyper{t,r},ESSseaten_hyper{t,r},'o','Color','k','MarkerFaceColor','k','MarkerSize',markersz,'LineWidth',lw);
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
    l=size(ESSseaten_hyper{t,i},2);
    if l>0
        p1=plot(radvals_hyper(i)*ones(l,1),ESSseaten_hyper{t,i},'o','Color','k','MarkerFaceColor','k','MarkerSize',markersz,'LineWidth',lw);
    end
    l=size(ESSsgettoeat_hyper{t,i},2);
    if l>0
        p2=plot(radvals_hyper(i)*ones(l,1),ESSsgettoeat_hyper{t,i},'s','Color','k','MarkerFaceColor','k','MarkerSize',markersz,'LineWidth',lw);
    end
    l=size(ESSsboth_hyper{t,i},2);
    if l>0
        p3=plot(radvals_hyper(i)*ones(l,1),ESSsboth_hyper{t,i},'d','Color','k','MarkerFaceColor','k','MarkerSize',markersz,'LineWidth',lw);
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

filename=strcat('/Users/eleanorbrush/Desktop/','test','.pdf');
print(filename,'-dpdf','-r300');