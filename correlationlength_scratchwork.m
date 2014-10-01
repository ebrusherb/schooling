positions=unifrnd(0,1,N,2);
d=squareform(pdist(positions));

%%
strategy=numneighbors*ones(N,1);
M=zeros(N);
for ind=1:N
    [~, order]=sort(d(ind,:));
    neighbors=order(2:strategy(ind)+1);
    M(ind,neighbors)=1/strategy(ind);
end
M(1:N+1:end)=-1; 

%%
corrvec=[];
distvec=[];
b=2;
radius=.1;
for i=1:N
    receiver=i;
    allreceivers=d(receiver,:)<=radius;
    beta=zeros(N,1);
    beta(allreceivers)=1;
    bbeta=b*beta;
    corrs=paircorrelations(M,bbeta);
    uppercorr=triu(corrs,1);
    upperdist=triu(d,1);
    corrvectoadd=uppercorr(~~uppercorr);
    distvectoadd=upperdist(~~upperdist);
%     [min(corrvectoadd) max(corrvectoadd)]
    corrvec=[corrvec; corrvectoadd; diag(corrs)];
    distvec=[distvec; distvectoadd; diag(d)];
%     corrvec=[corrvec;corrvectoadd];
%     distvec=[distvec;distvectoadd];
end
%%
figure
plot(distvec,corrvec,'ok')
%%
distbins=0:.01:1.5;
[~,w]=histc(distvec,distbins);
l=length(distbins);
avgcorr=zeros(1,l);
for i=1:l
    now=(w==i);
    avgcorr(i)=mean(corrvec(now));
end

i=sum(avgcorr>=0);
if i<l
%     corrlength=avgcorr(i)/(avgcorr(i)-avgcorr(i+1))*(distbins(i+1)-distbins(i))+distbins(i);
    corrlength=find(avgcorr<0,1,'first');
    corrlength=distbins(corrlength);
else corrlength=max(distbins);
end
corrlength

figure
plot(distbins,avgcorr/max(avgcorr))
hold on
plot(distbins,zeros(length(distbins),1))
plot(corrlength*ones(2,1),get(gca,'ylim'))
% hold off
