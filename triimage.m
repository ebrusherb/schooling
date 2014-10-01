function triimage(scaledvals,vals,coords,cornerlabels,mainlabel)

%%
textfontsz=20;
labfontsz=15;
figure
set(gcf,'Color','w')
set(gcf,'PaperUnits','inches')
v=get(gcf,'Position');
ratio=v(4)/v(3);
w=6.83;
h=w*ratio;

c=colormap;
spaced=linspace(0,1,size(c,1));
l=size(coords,1);

xoffset=.25;
yoffset=.2;

hold on
axis off
axis equal

for i=1:l
    p=coords(i,:);
    valnow=scaledvals(i);
    ind=sum(spaced<=valnow);
    if ind==size(c,1)
        valcol=c(end,:);
    else
        valcol=c(ind,:)+(valnow-spaced(ind))/(spaced(ind+1)-spaced(ind))*(c(ind+1,:)-c(ind,:));
    end
    fill([p(1),p(1)+cstep,p(1)+cstep,p(1),p(1)-cstep,p(1)-cstep],[p(2)-2/sqrt(3)*cstep,p(2)-1/sqrt(3)*cstep,p(2)+1/sqrt(3)*cstep,p(2)+2/sqrt(3)*cstep,p(2)+1/sqrt(3)*cstep,p(2)-1/sqrt(3)*cstep],valcol,'EdgeColor','none');
end

text(0-xoffset,-yoffset,cellstr(cornerlabels(1)),'Color','k','HorizontalAlignment','left','FontSize',textfontsz)
text(2+xoffset,-yoffset,cellstr(cornerlabels(2)),'Color','k','HorizontalAlignment','right','FontSize',textfontsz)
text(1,sqrt(3)+yoffset,cellstr(cornerlabels(3)),'HorizontalAlignment','center','Color','k','Clipping','off','FontSize',textfontsz)
% text(0,sqrt(3),'MUTUAL','HorizontalAlignment','center','Color','k','Clipping','off','FontSize',textfontsz)
% text(0,sqrt(3)-yoffset,'INFORMATION','HorizontalAlignment','center','Color','k','Clipping','off','FontSize',textfontsz)

    
mybar=colorbar;
tickvec=0:.2:1;
ticklabvec=tickvec*(max(vals)-min(vals))+min(vals);
ticklabvec=round(ticklabvec*10)/10;
set(mybar,'YTick',tickvec,'YTickLabel',ticklabvec,'FontSize',labfontsz)
yl=get(mybar,'YLabel');
v=get(yl,'Position');
set(yl,'String','Accuracy','FontSize',textfontsz,'Rotation',270,'Position',v+[3 0 0])

set(gcf,'PaperSize',[w h]);
set(gcf,'PaperPosition',[-w/10 -h/50 1.1*w 1.05*h]);
