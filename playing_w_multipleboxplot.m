data=horzcat(resident_scores,invader_scores);
Mlab={'Resident', 'Invader'};
colors=cols;

% Get sizes
M=size(data,2);
L=size(data,1);

% Calculate the positions of the boxes
positions=1:0.25:M*L*0.25+1+0.25*L;
positions(1:M+1:end)=[];

% Extract data and label it in the group correctly
x=[];
group=[];
toclear=[];
whichgroup=[];
for ii=1:L
    for jj=1:M
        aux=data{ii,jj};
        x=vertcat(x,aux(:));
        group=vertcat(group,ones(size(aux(:)))*jj+(ii-1)*M);
        if isempty(aux)
            toclear=[toclear (ii-1)*M+jj];
        else
            whichgroup=vertcat(whichgroup,jj);
        end
    end
end
shrunkpositions=positions;
shrunkpositions(toclear)=[];
% Plot it
boxplot(x,group, 'positions', shrunkpositions);

% Set the Xlabels
aux=reshape(positions,M,[]);
labelpos = sum(aux,1)./M;

set(gca,'xtick',labelpos)

    set(gca,'xticklabel',xlab);


% Get some colors

    cmap=colors;

color=repmat(cmap, 1, L);

% Apply colors
h = findobj(gca,'Tag','Box');
p = zeros(1,M);
for jj=1:length(h)
   p(whichgroup(length(h)-jj+1))=patch(get(h(jj),'XData'),get(h(jj),'YData'),color(1:3,jj)','FaceAlpha',color(4,jj));
end

legend(p,Mlab);

