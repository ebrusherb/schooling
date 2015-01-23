
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
%% lots of fitness and rho heatmaps

myseqmap=cbrewer('seq', 'YlOrRd',9);
myseqmap=myseqmap(end:-1:1,:);
% myseqmap=cbrewer('div','RdBu',9);
k=1;
ESSseaten=ESSeaten_total{k};
ESSsgettoeat=ESSgettoeat_total{k};
t=3;
radtoplot=[1 2 4:2:10 length(radvals)];
figure
for r=1:length(radtoplot)
    subplot(2,length(radtoplot),r)
    fitnesseaten=storefitnesseaten{t,radtoplot(r)};
    fitim=imagesc(transpose(fitnesseaten));
    set(gca,'ydir','normal')
    colormap(myseqmap)
%     colorbar
    caxis manual
    caxis([.8 1])
    title(strcat('Rad=',num2str(radvals(radtoplot(r)))));
    hold on
    plot(ESSseaten{t,radtoplot(r)},ESSseaten{t,radtoplot(r)},'or')
    subplot(2,length(radtoplot),r+length(radtoplot))
    rhoeaten=storerhoeaten{t,radtoplot(r)};
    rhoim=imagesc(transpose(rhoeaten));
    set(gca,'ydir','normal')
%     colormap(mydivmap)
%     colorbar
    caxis manual
    caxis([0 .1])
    hold on
    plot(ESSseaten{t,radtoplot(r)},ESSseaten{t,radtoplot(r)},'or')
end
%% one fitness and rho heatmap
t=3;
r=2;
figure
set(gcf,'Position',[200 200 1200 700])
subplot(1,2,1)
fitnesseaten=storefitnesseaten{t,radtoplot(r)};
    fitim=imagesc(transpose(fitnesseaten));
    set(gca,'ydir','normal')
    colormap(myseqmap)
    colorbar
    caxis manual
    caxis([.8 1])
    hold on
    plot(ESSseaten{t,radtoplot(r)},ESSseaten{t,radtoplot(r)},'or')
subplot(1,2,2)    
    rhoeaten=storerhoeaten{t,radtoplot(r)};
    rhoim=imagesc(transpose(rhoeaten));
    set(gca,'ydir','normal')
%     colormap(mydivmap)
    colorbar
    caxis manual
    caxis([.04 .06])
    hold on
    plot(ESSseaten{t,radtoplot(r)},ESSseaten{t,radtoplot(r)},'or')


%% ESS as a function of radius of signal
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

n=size(ESSseaten,2);
k=1;
t=3;

ESSseaten=ESSeaten_total{k};
ESSsgettoeat=ESSgettoeat_total{k};

hold on
for i=1:n
    l=size(ESSseaten{t,i},2);
    p1=plot(radvals(i)*ones(l,1),ESSseaten{t,i},'o','MarkerSize',markersz,'LineWidth',lw,'MarkerFaceColor','blue');
    l=size(ESSsgettoeat{t,i},2);
    if l>0
        p2=plot(radvals(i)*ones(l,1),ESSsgettoeat{t,i},'or','MarkerSize',markersz,'LineWidth',lw,'MarkerFaceColor','red');
    end
end
xlabel('Radius of signal','FontName',fontname,'FontSize',textfontsz)
ylabel('Number of neighbors','FontName',fontname,'FontSize',textfontsz)

leg=legend([p1 p2],{'Predation','Resources'},'FontName',fontname,'FontSize',labfontsz);
legend('boxoff')
v=get(leg,'Position');
set(leg,'Position',[.1 .4 v(3) v(4)])

set(gca,'FontName',fontname,'FontSize',labfontsz)


set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','ESSneighbors','.pdf');
% print(filename,'-dpdf','-r300');
%%  ESS as a function of radius for multiple times
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=.35*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
fontname='TimesNewRoman';
textfontsz=12;
labfontsz=10;
markersz=3;
lw=2;

n=size(ESSseaten,2);

for k=1:1
    ESSseaten=ESSeaten_total{k};
    ESSsgettoeat=ESSgettoeat_total{k};
for t=2:4
    subplot(1,3,t-1+(k-1)*4)
    hold on
    for i=1:n
        l=size(ESSseaten{t,i},2);
        p1=plot(radvals(i)*ones(l,1),ESSseaten{t,i},'o','MarkerSize',markersz,'LineWidth',lw,'MarkerFaceColor','blue');
        l=size(ESSsgettoeat{t,i},2);
        if l>0
            p2=plot(radvals(i)*ones(l,1),ESSsgettoeat{t,i},'or','MarkerSize',markersz,'LineWidth',lw,'MarkerFaceColor','red');
        end
    end
    if k==1 && t==2
            xlabel('Radius of signal','FontName',fontname,'FontSize',textfontsz)
            ylabel('Number of neighbors','FontName',fontname,'FontSize',textfontsz)
    end
    switch k
        case 1
            title(strcat('T=',num2str(Tvals(t))),'FontName',fontname,'FontSize',textfontsz)
    end
    if k==1 && t==2 
            leg=legend([p1 p2],{'Predation','Resources'},'FontName',fontname,'FontSize',labfontsz);
            legend('boxoff')
            v=get(leg,'Position');
            set(leg,'Position',[.1 .4 v(3) v(4)])
    end
    set(gca,'FontName',fontname,'FontSize',labfontsz)
end
end

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','ESSneighbors_multT','.pdf');
print(filename,'-dpdf','-r300');

%% greedyopt figures
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

nt=size(evolution_eaten,3);

subplot(2,3,1)
plot(1:nt,reshape(evolution_eaten(1,:,:),N,[]))
set(gca,'FontName',fontname,'FontSize',textfontsz)
xlabel('','FontName',fontname,'FontSize',textfontsz)
ylabel('Number of neighbors','FontName',fontname,'FontSize',textfontsz)
        
subplot(2,3,2)
plot(1:nt,1./groupconsensus_eaten)
set(gca,'FontName',fontname,'FontSize',textfontsz)
xlabel('','FontName',fontname,'FontSize',textfontsz)
ylabel('H2 Robustness','FontName',fontname,'FontSize',textfontsz)

subplot(2,3,3)
plot(1:nt,corrlengths_eaten)
set(gca,'FontName',fontname,'FontSize',textfontsz)
xlabel('','FontName',fontname,'FontSize',textfontsz)
ylabel('Correlation length','FontName',fontname,'FontSize',textfontsz)

subplot(2,3,4)
plot(1:nt,reshape(evolution_gettoeat(1,:,:),N,[]))
set(gca,'FontName',fontname,'FontSize',textfontsz)
subplot(2,3,5)
plot(1:nt,1./groupconsensus_gettoeat)
set(gca,'FontName',fontname,'FontSize',textfontsz)
subplot(2,3,6)
plot(1:nt,corrlengths_gettoeat)
set(gca,'FontName',fontname,'FontSize',textfontsz)

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=strcat('/Users/eleanorbrush/Desktop/','greedyoptneighbors','.pdf');
% print(filename,'-dpdf','-r300');


%%
figure
w=6.83;
h=.25*w;
set(gcf,'Units','inches');
set(gcf,'Position',[11.5 3 w h]);
subplot(1,2,1)
toplot=corrlengths;
% toplot=corrlengths./repmat([1 radvals(2:end)],Lh,1);
a=min(min(toplot))-1;
A=max(max(toplot))-1;
A=max(abs(a),abs(A))+1;
mycolors=cbrewer('div','RdBu',20);
imagesc(radvals,homogenstrats,toplot)
% caxis manual
% caxis([-A A])
% colormap(mycolors)
colorbar
title('Correlation Length')
xlabel('Radius of signal')
ylabel('Strategy')

subplot(1,2,2)
hold on
t=3;
n=size(ESSseaten,2);
for i=1:n
    l=length(ESSseaten{t,i});
    plot(corrlengths(5,i)*ones(l,1),ESSseaten{t,i},'o')
end
xlabel('Correlation length of highly connected group')
ylabel('ESS strategy')

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename='/Users/eleanorbrush/Desktop/corrlength_exploration.pdf';
% print(filename,'-dpdf','-r300')

%%
seqcols=colormap(jet(Nr));
figure
hold on

for i=1:Nr
    plot(ones(Lh,1)./(groupconsensus(:,i)./col(homogenstrats)),corrlengths(:,i),'-o','Color',seqcols(i,:),'LineWidth',2)
end
    