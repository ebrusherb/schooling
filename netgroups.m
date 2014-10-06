function groups = netgroups(M)

[vecm,valm]=eig(M);
valm=diag(valm);
s=intersect(find(sigfig(valm,14)==0),find(imag(valm)==0));
vecm=vecm(:,s);
ls=size(s,1);
Mbin=M;
Mbin(M==-1)=0;
Mbin(~~Mbin)=1;
g=sparse(Mbin);
[~,Ctotal]=graphconncomp(g,'Directed','true');
nc=max(Ctotal);
groupreps=zeros(1,nc);
for i=1:nc
    groupreps(i)=find(Ctotal==i,1,'first');
end
vecm=vecm(groupreps,:);
vecm=sigfig(vecm,14);

groups=cell(ls,1);
for i=1:ls
    groups{i}=find(abs(vecm(:,i))==max(abs(unique(vecm(:,i)))));
end

for i=1:(ls-1)
    for j=(i+1):ls
        if ~isempty(intersect(groups{i},groups{j}))
            groups{i}=union(groups{i},groups{j});
            groups{j}=[];
        end
    end
end

emptycells = cellfun(@isempty,groups);
groups(emptycells)=[];
ls=size(groups,1);
for i=1:nc
    count=0;
    for j=1:ls
        if ismember(i,groups{j})
            count=count+1;
        end
    end
    if count==0
        groups=[groups;[i]];
    end
end
