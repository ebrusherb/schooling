load /Users/eleanorbrush/Desktop/ESS_4regimes_hypogeom_nummoves=1000_numpermove=1.mat
ESSseaten_hyper=ESSseaten;
storefitnesseaten_hyper=storefitnesseaten;
storerhoeaten_hyper=storerhoeaten;
ESSsgettoeat_hyper=ESSsgettoeat;
storefitnessgettoeat_hyper=storefitnessgettoeat;
storerhogettoeat_hyper=storerhogettoeat;
ESSsboth_hyper=ESSsboth;
storefitnessboth_hyper=storefitnessboth;
storerhoboth_hyper=storerhoboth;
ESSsgenerous_hyper=ESSsgenerous;
storefitnessgenerous_hyper=storefitnessgenerous;
storerhogenerous_hyper=storerhogenerous;
radvals_hyper=radvals;

load /Users/eleanorbrush/Desktop/ESS_4regimes__nummoves=1000_numpermove=1.mat

%%
i=i+1;
fit1=storefitnessgettoeat_hyper{i};
fit2=storefitnessgettoeat{1,i};
rho1=storerhogettoeat_hyper{i};
rho2=storerhogettoeat{1,i};
Mf=max(max(max(fit1)),max(max(fit2)));
Mr=max(max(max(rho1)),max(max(rho2)));

subplot(2,2,1)

imagesc(strats,strats,transpose(fit1));set(gca,'ydir','normal');
hold on
plot(ESSsgettoeat_hyper{i},ESSsgettoeat_hyper{i},'ok','MarkerSize',8,'LineWidth',2)
caxis([0 Mf])
hold off

subplot(2,2,2)
% i=1;

imagesc(strats,strats,transpose(fit2));set(gca,'ydir','normal');
hold on
plot(ESSsgettoeat{1,i},ESSsgettoeat{1,i},'ok','MarkerSize',8,'LineWidth',2)
caxis([0 Mf])
hold off

subplot(2,2,3)

imagesc(strats,strats,transpose(rho1));set(gca,'ydir','normal');
caxis([0 Mr])

subplot(2,2,4)
% i=1;

imagesc(strats,strats,transpose(rho2));set(gca,'ydir','normal');
caxis([0 Mr])

%%
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