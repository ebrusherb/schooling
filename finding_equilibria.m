    %%
    %%how group performance depends on the strategy of the group and the
    %%tradeoff
    n=40;
    alphavals=0:.25:1;
    strats=2:2:n;
    groupshort=zeros(1,length(strats));
    grouplong=zeros(1,length(strats));
    groupperf=zeros(length(alphavals),length(strats));
    indshort=zeros(1,length(strats));
    indlong=zeros(1,length(strats));
    indperf=zeros(length(alphavals),length(strats));
    
    
    for s=1:length(strats)
            
        strat=strats(s);
        strategy=strat*ones(1,n);
        M=makenet(strategy);
    
        groupshort(s)=consensusspeed(M);
        grouplong(s)=(H2norm(M,'additive')/sqrt(n));
    
        v1=speeds(M,'both','mean');
        v2=vars(M,'additive','mean');
        indshort(s)=v1(1);
        indlong(s)=v2(1);
        
        for a=1:length(alphavals)
            alpha=alphavals(a);
            groupperf(a,s)=1/((1-alpha)*groupshort(s)+alpha*grouplong(s));
            indperf(a,s)=1/((1-alpha)*v1(1)+alpha*v2(1));
        end
    
    end
        
    
    
    HM=HeatMap(groupperf,'RowLabels',alphavals,'ColumnLabels',strats);
    P=HM.plot;
    colormap hot
    colorbar('Peer',P)
    caxis(P,[0 max(max(groupperf))])
    xlab=xlabel(P,'Strategy');
    ylab=ylabel(P,'Tradeoff');
    main=title(P,['Group Performance' ', N=' num2str(n)]);
    set([xlab ylab],'FontName','Arial','FontSize',12);
    set(main, 'FontName','Arial','FontSize',15);
    xpos=get(xlab,'Position');
    set(xlab,'Position',xpos+[0 .5 0]);
    figname=['/Users/eleanorbrush/Desktop/groupperfN', num2str(n) ,'.eps'];
    print(gcf,'-depsc',figname);
    
    grouprows=num2cell(groupperf,2);
    indrows=num2cell(indperf,2);
    
    normalizedgroupperf=cellfun(@(x) x/max(x),grouprows,'UniformOutput',false);
    normalizedgroupperf=cell2mat(normalizedgroupperf);
    
    HM=HeatMap(normalizedgroupperf,'RowLabels',alphavals,'ColumnLabels',strats);
    P=HM.plot;
    colormap hot
    colorbar('Peer',P)
    caxis(P,[0 max(max(normalizedgroupperf))])
    xlab=xlabel(P,'Strategy');
    ylab=ylabel(P,'Tradeoff');
    main=title(P,['Group Performance' ', N=' num2str(n)]);
    set([xlab ylab],'FontName','Arial','FontSize',12);
    set(main, 'FontName','Arial','FontSize',15);
    xpos=get(xlab,'Position');
    set(xlab,'Position',xpos+[0 1 0]);
    set(gca,'box','off')
    figname=['/Users/eleanorbrush/Desktop/normalizedgroupperfN', num2str(n) ,'.eps']; 
    print(gcf,'-depsc',figname);
    
    [optimalperf,~]=max(groupperf,[],2);
    groupoptstrats=cellfun(@(x) strats(x==max(x)),grouprows,'UniformOutput',false);
    indoptstrats=cellfun(@(x) strats(x==max(x)),indrows,'UniformOutput',false);
    if max(cellfun(@length,groupoptstrats))==1
        groupoptstrats=cell2mat(groupoptstrats);
    end
    if max(cellfun(@length,indoptstrats))==1
        indoptstrats=cell2mat(indoptstrats);
    end
%%
%if alpha=0, it's optimal to pay attention to as many neighbors as possible 
%if alpha=1, it's not optimal to pay attentiont to as few neighbors as
% possible-- how does the optimal number scale with the size of the group?
alpha=1;
nvals=20:20:40;
ln=length(nvals);

indoptstrats=cell(1,ln);
groupoptstrats=cell(1,ln);

for i=1:ln
    n=nvals(i);
    strats=2:2:n-1;
    groupshort=zeros(1,length(strats));
    grouplong=zeros(1,length(strats));
    groupperf=zeros(1,length(strats));
    indshort=zeros(1,length(strats));
    indlong=zeros(1,length(strats));
    indperf=zeros(1,length(strats));
    
    
    for s=1:length(strats)
            
        strat=strats(s);
        strategy=strat*ones(1,n);
        M=makenet(strategy);
    
        v1=speeds(M,'both','mean');
        v2=vars(M,'additive','mean');
        indshort(s)=v1(1);
        indlong(s)=v2(1);
        indperf(s)=1/((1-alpha)*v1(1)+alpha*v2(1));
        
        groupshort(s)=consensusspeed(M);
        grouplong(s)=(H2norm(M,'additive')/sqrt(n));
        groupperf(s)=1/((1-alpha)*groupshort(s)+alpha*grouplong(s));
        
    
    end
    
    indoptstrats{i}=strats(indperf==max(indperf));
    groupoptstrats{i}=strats(groupperf==max(groupperf));
    
end

if max(cellfun(@length,indoptstrats))==1
        indoptstrats=cell2mat(indoptstrats);
end
if max(cellfun(@length,groupoptstrats))==1
        groupoptstrats=cell2mat(groupoptstrats);
end

figure;
hold on
indplot=plot(nvals,indoptstrats);
groupplot=plot(nvals,groupoptstrats);
set(indplot,'LineWidth',2);
set(groupplot,'color','red','LineWidth',2);
xlab=xlabel('Size of Network');
ylab=ylabel('Optimal Strategy');
leg=legend([indplot groupplot],'External Signal','Consensus','location','best');
set([gca xlab ylab],'FontName','Arial','FontSize',12);
set(gcf,'PaperPositionMode','auto');
%figname=['/Users/eleanorbrush/Desktop/optstrats_alpha', num2str(alpha) ,'.eps']; 
%print(gcf,'-depsc',figname);

%%
% find ESS, NIS, possible dimorphisms, CSS
alphavals=[0];
la=length(alphavals);
Nvals=[40 80];
ln=length(Nvals);

ESSs=cell(ln,la);
NISs=cell(ln,la);
CSSs=cell(ln,la);
dimorphs=cell(ln,la);

for a=1:la
    for n=1:ln

        alpha=alphavals(a);
        N=Nvals(n);

    strats=2:2:(N-1);
    L=length(strats);
    fitness=zeros(L,L);
    rho=zeros(L,L);
    

        for i=1:L
            for j=1:L
            resident=strats(i);
            invader=strats(j);

            f=zeros(1,N-1);
            g=zeros(1,N-1);

            for k=1:(N-1)
                strategy=resident*ones(1,N);
                strategy(1:k)=invader;
                M=makenet(strategy);
                v1=speeds(M,'both','mean');
                v2=vars(M,'additive','mean');
                perf=1./((1-alpha)*v1+alpha*v2);
                f(k)=mean(perf(1:k));
                g(k)=mean(perf((k+1):N));
                if k==1
                    fitness(i,j)=f(k)/g(k);
                end
        
            end
            
            ratio=g./f;
            tosum=ones(1,N-1);
            for k=1:(N-1)
                tosum(k)=prod(ratio(1:k));
            end
            rho(i,j)=1/(1+sum(tosum));
            
            end
        
        end
        
       
            ESSs{n,a}=strats(ESSN(fitness,rho,N));
       
            NISs{n,a}=strats(NISN(fitness,rho,N));
        
            CSSs{n,a}=intersect(strats(ESSN(fitness,rho,N)),strats(NISN(fitness,rho,N)));
        
            dimorphs{n,a}=strats(dimorph(fitness));
    end  
end   

a=min(min(fitness));
b=max(max(fitness));
c=max(1-a,b-1);

%k=11;
n=size(fitness);
n=n(2);
F=fitness;
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

% figname=['/Users/eleanorbrush/Desktop/fitness_N',num2str(N),'_alpha', num2str(alpha) ,'.eps']; 
% print(gcf,'-depsc',figname);
%%  find the fraction of invaders that can invade
alpha=.25;
n=40;

range=1:n;
strats=2:2:(n-1);
L=length(strats);

rho=zeros(L,L);
fraction=zeros(L,L);    
        for i=1:L
            for j=1:L
            resident=strats(i);
            invader=strats(j);

            f=zeros(1,n-1);
            g=zeros(1,n-1);

            for k=1:(n-1)
                strategy=resident*ones(1,n);
                strategy(1:k)=invader;
                M=makenet(strategy);
                v1=speeds(M,'both','mean');
                v2=vars(M,'additive','mean');
                perf=1./((1-alpha)*v1+alpha*v2);
                f(k)=mean(perf(1:k));
                g(k)=mean(perf((k+1):n));
        
            end

        ratio=g./f;
        tosum=ones(1,n-1);
        for k=1:(n-1)
            tosum(k)=prod(ratio(1:k));
        end
        rho(i,j)=1/(1+sum(tosum));
        
        probup=((1:(n-1)).*f)./((1:(n-1)).*f+((n-1):-1:1).*g).*(((n-1):-1:1))/n;
            probdown=(((n-1):-1:1).*g)./((1:(n-1)).*f+((n-1):-1:1).*g).*(1:(n-1))/n;
            if i~=j
                if sum((probup-probdown)<0)==0
                    fraction(i,j)=n;
                else
                fraction(i,j)=min(range((probup-probdown)<0))-1;
                if fraction(i,j)==n-1
                    fraction(i,j)=n;
                end
                if fraction(i,j)==-1
                    fraction(i,j)=0;
                end
                end
            end
        
        
            end
        end
 HM=HeatMap(transpose(rho-1/n),'RowLabels',strats,'ColumnLabels',strats);
    addXLabel(HM,'Resident Strategy');
    addYLabel(HM,'Invader Strategy');
    P=HM.plot;
    %colormap hot
    colorbar('Peer',P)
    caxis(P,[-.1 .1])       
        
HM=HeatMap(transpose(fraction),'RowLabels',strats,'ColumnLabels',strats);
    P=HM.plot;
    colormap jet
    colorbar('Peer',P)
    caxis(P,[0 max(max(fraction))])
    xlab=xlabel(P,'Strategy');
    ylab=ylabel(P,'Tradeoff');
    set([xlab ylab],'FontName','Arial','FontSize',12);
    xpos=get(xlab,'Position');
    set(xlab,'Position',xpos+[0 1 0]);
    set(gca,'box','off')
figname=['/Users/eleanorbrush/Desktop/numinvaders_N',num2str(n),'_alpha', num2str(alpha) ,'.eps']; 
print(gcf,'-depsc',figname);

    %%
    %how group performance is improved or decreased when invaders invade as much as they can

    
    n=40;
    alpha=.25;
    strats=2:2:(n-1);
    L=length(strats);
    
    groupshort=zeros(L,L,n+1);
    grouplong=zeros(L,L,n+1);
    residentperf=zeros(L,L,n+1);
    invaderperf=zeros(L,L,n+1);
    groupperf=zeros(L,L,n+1);
    groupimprovement=zeros(L,L);
    
    for I=1:length(strats)
        for J=1:length(strats)
    
    resident=strats(I);
    invader=strats(J);
    
    strategy=resident*ones(1,n);
    M=makenet(strategy);
    v1=speeds(M,'both','mean');
    v2=vars(M,'additive','mean');
    perf=1./((1-alpha)*v1+alpha*v2);
    residentperf(I,J,1)=mean(perf);
    grouplong(I,J,1)=(H2norm(M,'additive')/sqrt(n));
    groupshort(I,J,1)=consensusspeed(M);
    groupperf(I,J,1)=1/((1-alpha)*groupshort(I,J,1)+alpha*grouplong(I,J,1));
    
    for k=1:(n-1)
        strategy=resident*ones(1,n);
        strategy(1:k)=invader;
        M=makenet(strategy);
        v1=speeds(M,'both','mean');
        v2=vars(M,'additive','mean');
        perf=1./((1-alpha)*v1+alpha*v2);
        residentperf(I,J,k+1)=mean(perf((k+1):n));
        invaderperf(I,J,k+1)=mean(perf(1:k));
        grouplong(I,J,k+1)=(H2norm(M,'additive')/sqrt(n));
        groupshort(I,J,k+1)=consensusspeed(M);
        groupperf(I,J,k+1)=1/((1-alpha)*groupshort(I,J,k+1)+alpha*grouplong(I,J,k+1));
    end
    
    strategy=invader*ones(1,n);
    M=makenet(strategy);
    v1=speeds(M,'both','mean');
    v2=vars(M,'additive','mean');
    perf=1./((1-alpha)*v1+alpha*v2);
    invaderperf(I,J,n+1)=mean(perf);
    grouplong(I,J,n+1)=(H2norm(M,'additive')/sqrt(n));
    groupshort(I,J,n+1)=consensusspeed(M);
    groupperf(I,J,n+1)=1/((1-alpha)*groupshort(I,J,n+1)+alpha*grouplong(I,J,n+1));
%figure   
%v=invaderperf(1:N);
%plot(0:(N-1),v)
%hold all
%v=residentperf(1:N);
%plot(0:(N-1),v)
%v=groupperf(1:N);plot(0:(N-1),v)
%plot(fraction(I,J),-4:.1:3,'ro')
%     groupimprovement(I,J)=groupperf(1)-groupperf(fraction(I,J)+1);
        end
    end
  
HM=HeatMap(transpose(groupimprovement),'RowLabels',strats,'ColumnLabels',strats);
    addXLabel(HM,'Resident Strategy');
    addYLabel(HM,'Invader Strategy');
    addTitle(HM,['Improvement in Group Performance' ', alpha =' num2str(alpha)]);
    P=HM.plot;
    %colormap hot
    colorbar('Peer',P)
    m=max(abs(min(min(groupimprovement))),abs(max(max(groupimprovement))));
    caxis(P,[-m m])
    

    
 