function [distbins,avgcorr,corrlength] = correlationlength_mat_single(M,d,b,radius,i)
% if nargin==5 || strcmp(path,'path')
%     path=1;
% else path=0;
% end
% 
% if ~path
N=size(M,2);
corrvec=[];
distvec=[];

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
corrvec=[corrvec; corrvectoadd; diag(corrs)];
distvec=[distvec; distvectoadd; diag(d)];

distbins=0:.02:1.5;
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
else corrlength=max(distvec);
end
% end
% 
% if path
% M(M==-1)=0;
% M(~~M)=1;
% g=sparse(M);
% paths=graphallshortestpaths(g);
% 
% 
% corrvec=[];
% pathvec=[];
% for i=1:N
%     receiver=i;
%     allreceivers=d(receiver,:)<=radius;
%     beta=zeros(N,1);
%     beta(allreceivers)=1;
%     bbeta=b*beta;
%     corrs=paircorrelations(M,bbeta);
%     uppercorr=triu(corrs,1);
%     upperpath=triu(paths,1);
%     lowerpath=triu(transpose(paths),1);
%     corrvectoadd=uppercorr(~~uppercorr);
%     pathvectoadd=[upperpath(~~upperpath) lowerpath(~~lowerpath)];
%     pathvectoadd=min(pathvectoadd')';
%     corrvec=[corrvec; corrvectoadd; diag(corrs)];
%     pathvec=[pathvec; pathvectoadd; diag(paths)];
% end
% pathbins=0:1:max(max(paths));
% 
% [~,w]=histc(pathvec,pathbins);
% l=length(pathbins);
% avgcorr=zeros(1,l);
% for i=1:l
%     now=(w==i);
%     avgcorr(i)=mean(corrvec(now));
% end
% i=sum(avgcorr>=0);
% if i<l
%     corrlength=avgcorr(i)/(avgcorr(i)-avgcorr(i+1))*(distbins(i+1)-distbins(i))+distbins(i);
% else corrlength=max(pathvec);
% end
% end
% 
% end