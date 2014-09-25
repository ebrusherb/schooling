%%
axislw=1.25;
lw=3;
labfontsz=10 ;
textfontsz=12;
alpha=.7;
ticklength=[.02,.0250];
markersz=5;
markercol='white';
fontname='Times New Roman';
%%
% find ESS, NIS, possible dimorphisms, CSS
numsigs_permove=1;
nummoves=1000;
N=20;
b=1;

Tvals=[1];
NT=length(Tvals);
% radvals=[.1 .5 1];
radvals=.25;
Nr=length(radvals);

ESSseaten=cell(NT,Nr);
NISseaten=cell(NT,Nr);
CSSseaten=cell(NT,Nr);

ESSsgettoeat=cell(NT,Nr);
NISsgettoeat=cell(NT,Nr);
CSSsgettoeat=cell(NT,Nr);

strats=2:2:(N-1);
L=length(strats);

for i=1:NT
    for j=1:Nr
        T=Tvals(i);
        radius=radvals(j);
        fitnesseaten=zeros(L,L);
        rhoeaten=zeros(L,L);
        
        fitnessgettoeat=zeros(L,L);
        rhogettoeat=zeros(L,L);


        for u=1:L
            for v=1:L
            resident=strats(u);
            invader=strats(v);

            featen=zeros(1,N-1);
            geaten=zeros(1,N-1);
            
            fgettoeat=zeros(1,N-1);
            ggettoeat=zeros(1,N-1);

            for k=1:(N-1)
                strategy=resident*ones(1,N);
                strategy(1:k)=invader;
                [meanscore scorevar probeaten probgettoeat]=signalingevents(strategy,numsigs_permove,nummoves,radius,b,T);
                perfeaten=1-probeaten;
                featen(k)=mean(perfeaten(1:k));
                geaten(k)=mean(perfeaten((k+1):N));
                if k==1
                    fitnesseaten(u,v)=featen(k)/geaten(k);
                end
                
                perfgettoeat=probgettoeat;
                fgettoeat(k)=mean(perfgettoeat(1:k));
                ggettoeat(k)=mean(perfgettoeat((k+1):N));
                if k==1
                    fitnessgettoeat(u,v)=fgettoeat(k)/ggettoeat(k);
                end 
            end

            ratio=geaten./featen;
            tosum=ones(1,N-1);
            for k=1:(N-1)
                tosum(k)=prod(ratio(1:k));
            end
            rhoeaten(u,v)=1/(1+sum(tosum));
            
            ratio=ggettoeat./fgettoeat;
            tosum=ones(1,N-1);
            for k=1:(N-1)
                tosum(k)=prod(ratio(1:k));
            end
            rhogettoeat(u,v)=1/(1+sum(tosum));

            end

        end

        eqstratseaten=eq_strats(N,fitnesseaten,rhoeaten);
        ESSseaten{i,j}=eqstratseaten{2};
        NISseaten{i,j}=eqstratseaten{3};
        CSSseaten{i,j}=eqstratseaten{4};
        
        eqstratsgettoeat=eq_strats(N,fitnessgettoeat,rhogettoeat);
        ESSsgettoeat{i,j}=eqstratsgettoeat{2};
        NISsgettoeat{i,j}=eqstratsgettoeat{3};
        CSSsgettoeat{i,j}=eqstratsgettoeat{4};
    end
end

%%
a=min(min(fitnesseaten));
b=max(max(fitnesseaten));
c=max(1-a,b-1);

%%
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=.75*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
% set(gca,'FontSize',labfontsz)

n=size(fitnesseaten);
n=n(2);
F=fitnesseaten;
%F(logical(eye(size(F))))=a-2*(b-a)/k;
I=imagesc(transpose(F(:,n:-1:1)));
%cols=colormap(jet(k));
%cols=[[0 0 0];cols];
%set(gcf,'colormap',cols)
bar=colorbar;
set(gca,'XTick',1:2:length(strats),'XTickLabel',strats(1:2:length(strats)))
set(gca,'YTick',1:2:length(strats),'YTickLabel',strats(length(strats):-2:1))
xlab=xlabel('Resident Strategy');
ylab=ylabel('Invader Strategy');
set(get(bar,'ylabel'),'string','Fitness','FontName','Arial','FontSize',12);
set([xlab ylab],'FontName','Arial','FontSize',12);

% hold on
% 
% v=NISN(fitness,rho,N);
% for i=1:length(v)
%     plot(0:(length(strats)+1),(length(strats)-v(i)+1)*ones(length(strats)+2,1),'-','color','white')
% end

% filename=['/Users/eleanorbrush/Desktop/fitness_N',num2str(N),'_rad', num2str(radius) ,'.pdf']; 
% print(filename,'-dpdf','-r300');

%%
textfontsz=20;
labfontsz=20;
figure
set(gcf,'Color','w')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=.75*w;
set(gcf,'Units','inches');
set(gcf,'Position',[.5 1 w h]);
set(gca,'FontSize',labfontsz)


hold on
i=1;
    plot(radvals(i)*ones(length(ESSseaten{1,i}),1),strats(ESSseaten{1,i}),'or','LineWidth',lw,'MarkerFaceColor','red')
i=2;
    plot(radvals(i),strats(ESSseaten{1,i}(1)),'or','LineWidth',lw,'MarkerFaceColor','red')
i=3;
    plot(radvals(i)*ones(length(ESSseaten{1,i}),1),strats(ESSseaten{1,i}),'or','LineWidth',lw,'MarkerFaceColor','red')
    
xlabel('Radius of signal','FontSize',textfontsz);
ylabel('ESS # of neighbors','FontSize',textfontsz);

set(gca,'ylim',[0 20])
set(gca,'xlim',[0 1])

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[0 0 w h]);

filename=['/Users/eleanorbrush/Desktop/fitness_N',num2str(N),'_rad', num2str(radius) ,'.pdf']; 
print(filename,'-dpdf','-r300');


%%
% %%  find the fraction of invaders that can invade
% alpha=.25;
% n=40;
% 
% range=1:n;
% strats=2:2:(n-1);
% L=length(strats);
% 
% rho=zeros(L,L);
% fraction=zeros(L,L);    
%         for i=1:L
%             for j=1:L
%             resident=strats(i);
%             invader=strats(j);
% 
%             f=zeros(1,n-1);
%             g=zeros(1,n-1);
% 
%             for k=1:(n-1)
%                 strategy=resident*ones(1,n);
%                 strategy(1:k)=invader;
%                 M=makenet(strategy);
%                 v1=speeds(M,'both','mean');
%                 v2=vars(M,'additive','mean');
%                 perf=1./((1-alpha)*v1+alpha*v2);
%                 f(k)=mean(perf(1:k));
%                 g(k)=mean(perf((k+1):n));
%         
%             end
% 
%         ratio=g./f;
%         tosum=ones(1,n-1);
%         for k=1:(n-1)
%             tosum(k)=prod(ratio(1:k));
%         end
%         rho(i,j)=1/(1+sum(tosum));
%         
%         probup=((1:(n-1)).*f)./((1:(n-1)).*f+((n-1):-1:1).*g).*(((n-1):-1:1))/n;
%             probdown=(((n-1):-1:1).*g)./((1:(n-1)).*f+((n-1):-1:1).*g).*(1:(n-1))/n;
%             if i~=j
%                 if sum((probup-probdown)<0)==0
%                     fraction(i,j)=n;
%                 else
%                 fraction(i,j)=min(range((probup-probdown)<0))-1;
%                 if fraction(i,j)==n-1
%                     fraction(i,j)=n;
%                 end
%                 if fraction(i,j)==-1
%                     fraction(i,j)=0;
%                 end
%                 end
%             end
%         
%         
%             end
%         end
%  HM=HeatMap(transpose(rho-1/n),'RowLabels',strats,'ColumnLabels',strats);
%     addXLabel(HM,'Resident Strategy');
%     addYLabel(HM,'Invader Strategy');
%     P=HM.plot;
%     %colormap hot
%     colorbar('Peer',P)
%     caxis(P,[-.1 .1])       
%         
% HM=HeatMap(transpose(fraction),'RowLabels',strats,'ColumnLabels',strats);
%     P=HM.plot;
%     colormap jet
%     colorbar('Peer',P)
%     caxis(P,[0 max(max(fraction))])
%     xlab=xlabel(P,'Strategy');
%     ylab=ylabel(P,'Tradeoff');
%     set([xlab ylab],'FontName','Arial','FontSize',12);
%     xpos=get(xlab,'Position');
%     set(xlab,'Position',xpos+[0 1 0]);
%     set(gca,'box','off')
% figname=['/Users/eleanorbrush/Desktop/numinvaders_N',num2str(n),'_alpha', num2str(alpha) ,'.eps']; 
% print(gcf,'-depsc',figname);
% 
    %%
%     %how group performance is improved or decreased when invaders invade as much as they can
% 
%     
%     n=40;
%     alpha=.25;
%     strats=2:2:(n-1);
%     L=length(strats);
%     
%     groupshort=zeros(L,L,n+1);
%     grouplong=zeros(L,L,n+1);
%     residentperf=zeros(L,L,n+1);
%     invaderperf=zeros(L,L,n+1);
%     groupperf=zeros(L,L,n+1);
%     groupimprovement=zeros(L,L);
%     
%     for I=1:length(strats)
%         for J=1:length(strats)
%     
%     resident=strats(I);
%     invader=strats(J);
%     
%     strategy=resident*ones(1,n);
%     M=makenet(strategy);
%     v1=speeds(M,'both','mean');
%     v2=vars(M,'additive','mean');
%     perf=1./((1-alpha)*v1+alpha*v2);
%     residentperf(I,J,1)=mean(perf);
%     grouplong(I,J,1)=(H2norm(M,'additive')/sqrt(n));
%     groupshort(I,J,1)=consensusspeed(M);
%     groupperf(I,J,1)=1/((1-alpha)*groupshort(I,J,1)+alpha*grouplong(I,J,1));
%     
%     for k=1:(n-1)
%         strategy=resident*ones(1,n);
%         strategy(1:k)=invader;
%         M=makenet(strategy);
%         v1=speeds(M,'both','mean');
%         v2=vars(M,'additive','mean');
%         perf=1./((1-alpha)*v1+alpha*v2);
%         residentperf(I,J,k+1)=mean(perf((k+1):n));
%         invaderperf(I,J,k+1)=mean(perf(1:k));
%         grouplong(I,J,k+1)=(H2norm(M,'additive')/sqrt(n));
%         groupshort(I,J,k+1)=consensusspeed(M);
%         groupperf(I,J,k+1)=1/((1-alpha)*groupshort(I,J,k+1)+alpha*grouplong(I,J,k+1));
%     end
%     
%     strategy=invader*ones(1,n);
%     M=makenet(strategy);
%     v1=speeds(M,'both','mean');
%     v2=vars(M,'additive','mean');
%     perf=1./((1-alpha)*v1+alpha*v2);
%     invaderperf(I,J,n+1)=mean(perf);
%     grouplong(I,J,n+1)=(H2norm(M,'additive')/sqrt(n));
%     groupshort(I,J,n+1)=consensusspeed(M);
%     groupperf(I,J,n+1)=1/((1-alpha)*groupshort(I,J,n+1)+alpha*grouplong(I,J,n+1));
% %figure   
% %v=invaderperf(1:N);
% %plot(0:(N-1),v)
% %hold all
% %v=residentperf(1:N);
% %plot(0:(N-1),v)
% %v=groupperf(1:N);plot(0:(N-1),v)
% %plot(fraction(I,J),-4:.1:3,'ro')
% %     groupimprovement(I,J)=groupperf(1)-groupperf(fraction(I,J)+1);
%         end
%     end
%   
% HM=HeatMap(transpose(groupimprovement),'RowLabels',strats,'ColumnLabels',strats);
%     addXLabel(HM,'Resident Strategy');
%     addYLabel(HM,'Invader Strategy');
%     addTitle(HM,['Improvement in Group Performance' ', alpha =' num2str(alpha)]);
%     P=HM.plot;
%     %colormap hot
%     colorbar('Peer',P)
%     m=max(abs(min(min(groupimprovement))),abs(max(max(groupimprovement))));
%     caxis(P,[-m m])
%     
%%
% N=40;
% beta=-.1;
% strategy=32*ones(1,N);
% 
% strats=2:2:(N-1);
% l=length(strats);
% v=zeros(1,l);
% perf=zeros(N,l);
% 
% for i=1:l
%     strategy(1)=strats(i);
%     c=covariance(strategy,beta);
%     perf(:,i)=c{3};
% end
%  