function [distvec,corrvec,corrlength] = correlationlength_mat_single_v2(M,d,b,radius,i)
% if nargin==5 || strcmp(path,'path')
%     path=1;
% else path=0;
% end
% 
% if ~path
N=size(M,2);

receiver=i;
allreceivers=d(receiver,:)<=radius;
beta=zeros(N,1);
beta(allreceivers)=1;
bbeta=b*beta;
corrs=paircorrelations(M,bbeta);

corrvec=corrs(i,:);
distvec=d(i,:);
[~,o]=sort(distvec);
distvec=distvec(o);
corrvec=corrvec(o);

l=length(distvec);

i=sum(corrvec>=0);
if i<l
    corrlength=find(corrvec<0,1,'first');
    corrlength=distvec(corrlength);
else corrlength=max(distvec);
end
