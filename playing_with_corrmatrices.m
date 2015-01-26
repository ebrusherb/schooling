m=[3 1;4 3];
% theta=pi/3;
% m=[cos(theta) -sin(theta); sin(theta) cos(theta)];
minv=inv(m);

minv=[4 .01;6 5];
m=inv(minv);
[vecs,vals]=eig(minv);
vals=diag(vals);
[1/vals(1) 1/vals(2)]
v1=vecs(:,1);
v2=vecs(:,2);
Linv=vecs;
L=inv(Linv);
v=[4;5];
c=L*v;

xvals=-5:.1:5;
yvals=-5:.1:5;
nx=length(xvals);
ny=length(yvals);

tocontour=zeros(nx,ny);

for i=1:nx
    for j=1:ny
        x=xvals(i);
        y=yvals(j);
        vnow=[x; y];
        tocontour(i,j)=transpose(vnow)*minv*vnow;
    end
end

imagesc(xvals,yvals,tocontour)
colorbar
axis square
set(gca,'ydir','normal')
% contour(repmat(xvals,ny,1),repmat(col(yvals),1,nx),tocontour)
hold on
plot([0 v1(1)],[0 v1(2)],'r','LineWidth',2)
plot([0 v2(1)],[0 v2(2)],'r','LineWidth',2)
hold off