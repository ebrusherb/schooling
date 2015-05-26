N=40;
m=4;
resident=N-1;

probeaten=zeros(N-1,2);
toadd_keep=zeros(N-1,m+1);

for j=1:(N-2)
    invader=j;
    k=1;

    %             strategy=resident*ones(1,N);
    %             strategy(N+1-(1:k))=invader;
    %             [probeaten, probgettoeat, both, generous]=signalingevents_hypergeo_parallel(strategy,numsigs_permove,nummoves,m);

    toadd=zeros(m+1,1);
    for l=0:m
       toadd(l+1)=hygepdf(l,N-1,m,invader)*power(1-hygecdf(l,N-1,m,invader),k-1)*power(1-hygecdf(floor(l/invader*resident),N-1,m,resident),N-k);
    end
    toadd_keep(j,:)=toadd;
    probleast=sum(toadd);
%     probeaten=zeros(1,2);
    probeaten(j,1)=probleast;
    probeaten(j,2)=(1-probleast*k)/(N-k);
end

%%
N=20;
k=10;
resident=N-3;
invader=N-5;

toadd1=zeros(m+1,1);
for l=0:m
   toadd1(l+1)=hygepdf(l,N-1,m,invader)*power(1-hygecdf(l,N-1,m,invader),k-1)*power(1-hygecdf(floor(l/invader*resident),N-1,m,resident),N-k);
end
p1=sum(toadd1);

toadd2=zeros(m+1,1);
for l=0:m
   toadd2(l+1)=hygepdf(l,N-1,m,resident)*power(1-hygecdf(l,N-1,m,invader),N-k-1)*power(1-hygecdf(floor(l/resident*invader),N-1,m,invader),k);
end
p2=sum(toadd2);