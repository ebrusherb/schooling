function [probeaten, probgettoeat, both, probgettoeat_generous]=signalingevents_hypergeo_parallel(strategy,numsigs_permove,nummoves,m)
N=max(size(strategy));

numsigs_tot=numsigs_permove*nummoves;


scores=zeros(N,numsigs_permove,nummoves);

parfor i=1:nummoves
% for i=1:nummoves
    
    for j=1:numsigs_permove
        draw=hygernd(N,m,strategy);
        freq=draw./strategy;
        
        scores(:,j,i)=freq;
    end
end

scores=reshape(scores,N,[]);

[minvals,~]=min(scores);
minmat=repmat(minvals,N,1);
[rows,cols]=find(abs(scores-minmat)<0.00001);
minscorer=zeros(N,size(scores,2));

for i=1:size(scores,2)
    look=find(cols==i);
    minscorer(rows(look),i)=1/size(look,1)/numsigs_tot;
%     minscorer(rows(cols==i),i)=1/numsigs_tot;
end
probeaten=sum(minscorer,2);

[maxvals,~]=max(scores);
maxmat=repmat(maxvals,N,1);
[rows,cols]=find(abs(scores-maxmat)<0.00001);
maxscorer=zeros(N,size(scores,2));

for i=1:size(scores,2)
    look=find(cols==i);
    maxscorer(rows(look),i)=1/size(look,1)/numsigs_tot;
%     maxscorer(rows(cols==i),i)=1/numsigs_tot;
end
probgettoeat=sum(maxscorer,2);
% end

scores_first=scores(:,1:nummoves/2);
scores_next=scores(:,(nummoves/2+1):end);

[minvals_first,~]=min(scores_first);
minmat_first=repmat(minvals_first,N,1);
[rows_first,cols_first]=find(abs(scores_first-minmat_first)<0.00001);
minscorer_first=zeros(N,size(scores_first,2));

[maxvals_next,~]=max(scores_next);
maxmat_next=repmat(maxvals_next,N,1);
[rows_next,cols_next]=find(abs(scores_next-maxmat_next)<0.00001);
maxscorer_next=zeros(N,size(scores_next,2));

for i=1:size(scores_first,2)
    look=find(cols_first==i);
    minscorer_first(rows_first(look),i)=1/size(look,1)/(numsigs_tot/2);
end
for i=1:size(scores_next,2)
    look=find(cols_next==i);
    maxscorer_next(rows_next(look),i)=1/size(look,1)/(numsigs_tot/2);
end

probeaten_first=sum(minscorer_first,2);
probgettoeat_next=sum(maxscorer_next,2);
both=(1-probeaten_first).*probgettoeat_next;

cutoff=3;
maxscorer_generous=zeros(N,size(scores,2));
for i=1:size(scores,2)
    [x,~,z]=unique(scores(:,i));
    [~,o]=sort(x);
    cutoff_now=min(cutoff,length(x)-1);
%     maxscorer_generous(o(z)>=o(end-cutoff_now),i)=1/numsigs_tot;
    choose=zeros(N,1);
    count=0;
    for l=1:length(o)
        w=find(z==o(l));
        if length(w)>1
        choose(count+(1:length(w)))=randsample(w,length(w),'false');
        else choose(count+1)=w;
        end
        count=count+length(w);
    end
    maxscorer_generous(choose(end:-1:end-cutoff),i)=1/numsigs_tot;
end
probgettoeat_generous=sum(maxscorer_generous,2);