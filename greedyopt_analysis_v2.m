% load /Users/eleanorbrush/Desktop/greedyopt_rad=0.1_T=1_nummoves=1000_numpermove=1_timesteps=200.mat


%%

% evolution_eaten=zeros(its,N,timesteps)
its=size(evolution_eaten,1);
timesteps=size(evolution_eaten,3);
hold on

plot(1:its,var(evolution_eaten(:,:,end),[],2))
plot(1:its,var(evolution_gettoeat(:,:,end),[],2),'red')
legend({'Predation','Resources'})

figure
hist(col(evolution_eaten(:,:,end)))

%%
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
plot(1:timesteps,reshape(evolution_eaten(i,:,:),N,[]),'LineWidth',thinlw)
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
plot(1:timesteps,reshape(evolution_gettoeat(1,:,:),N,[]),'LineWidth',thinlw)
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
plot(1:timesteps,reshape(evolution_both(i,:,:),N,[]),'LineWidth',thinlw)
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
plot(1:timesteps,reshape(evolution_generous(i,:,:),N,[]),'LineWidth',thinlw)
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