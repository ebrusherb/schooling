
ESSeaten_total=cell(2,1);
ESSgettoeat_total=cell(2,1);

load /Users/eleanorbrush/Desktop/ESS_nummoves=1000_numpermove=1.mat

ESSeaten_total{1}=ESSseaten;
ESSgettoeat_total{1}=ESSsgettoeat;

load /Users/eleanorbrush/Desktop/ESS_nummoves=1_numpermove=1000.mat

ESSeaten_total{2}=ESSseaten;
ESSgettoeat_total{2}=ESSsgettoeat;
%%
k=1;
ESSseaten=ESSeaten_total{k};
ESSsgettoeat=ESSgettoeat_total{k};
mycolormap=cbrewer('seq', 'YlGnBu',9);
figure
imagesc(transpose(fitnesseaten))
set(gca,'ydir','normal')
colormap(mycolormap)
colorbar
ylabel('Invader')
xlabel('Resident')

%%
numsigs_permove=1;
nummoves=1000;

N=20;
b=1;
radius=radvals(end);
T=Tvals(4);

strategy=[19,16*ones(1,N-1)];
[probeaten, ~]=signalingevents(strategy,numsigs_permove,nummoves,radius,b,T);
perfeaten=1-probeaten
%%
figure
n=size(ESSseaten,2);

for k=1:2
    ESSseaten=ESSeaten_total{k};
    ESSsgettoeat=ESSgettoeat_total{k};
for t=1:4
    subplot(2,4,t+(k-1)*4)
    hold on
    for i=1:n
        l=size(ESSseaten{t,i},2);
        plot(radvals(i)*ones(l,1),ESSseaten{t,i},'o');
        l=size(ESSsgettoeat{t,i},2);
        plot(radvals(i)*ones(l,1),ESSsgettoeat{t,i},'or');
    end
    switch t+(k-1)*4
        case 1
            xlabel('Radius of signal')
            ylabel('Number of neighbors')
    end
    switch k
        case 1
            title(strcat('T=',num2str(Tvals(t))))
    end
end
end