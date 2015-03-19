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
h=.5*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
fontname='TimesNewRoman';
textfontsz=12;
markersz=3;
lw=2;

timesteps=size(evolution_eaten,3);
M=max([max(maxt_eaten),max(maxt_gettoeat)]);
M=min(timesteps,M+2);

subplot(2,3,1)
plot(1:timesteps,reshape(evolution_eaten(1,:,:),N,[]))
set(gca,'xlim',[0 M])
set(gca,'FontName',fontname,'FontSize',textfontsz)
xlabel('','FontName',fontname,'FontSize',textfontsz)
ylabel('Number of neighbors','FontName',fontname,'FontSize',textfontsz)
        
subplot(2,3,2)
plot(1:timesteps,1./groupconsensus_eaten)
set(gca,'xlim',[0 M])
set(gca,'FontName',fontname,'FontSize',textfontsz)
xlabel('','FontName',fontname,'FontSize',textfontsz)
ylabel('H2 Robustness','FontName',fontname,'FontSize',textfontsz)

subplot(2,3,3)
plot(1:timesteps,corrlengths_eaten)
set(gca,'xlim',[0 M])
set(gca,'FontName',fontname,'FontSize',textfontsz)
xlabel('','FontName',fontname,'FontSize',textfontsz)
ylabel('Correlation length','FontName',fontname,'FontSize',textfontsz)

subplot(2,3,4)
plot(1:timesteps,reshape(evolution_gettoeat(1,:,:),N,[]))
set(gca,'xlim',[0 M])
set(gca,'FontName',fontname,'FontSize',textfontsz)
subplot(2,3,5)
plot(1:timesteps,1./groupconsensus_gettoeat)
set(gca,'xlim',[0 M])
set(gca,'FontName',fontname,'FontSize',textfontsz)
subplot(2,3,6)
plot(1:timesteps,corrlengths_gettoeat)
set(gca,'xlim',[0 M])
set(gca,'FontName',fontname,'FontSize',textfontsz)

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','greedyoptneighbors','.pdf');
% print(filename,'-dpdf','-r300');