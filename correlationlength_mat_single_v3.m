function [distvec,corrvec,corrlength,crosscorrs] = correlationlength_mat_single_v3(M,d,b,radius,receiver)

N=size(M,2);

allreceivers=d(receiver,:)<=radius;
numreceivers=sum(allreceivers);
beta=zeros(N,1);
beta(allreceivers)=1;
bbeta=b*beta;
corrs=paircorrelations(M,bbeta);

corrvec=corrs(receiver,:);
distvec=d(receiver,:);
[~,o]=sort(distvec);
distvec=distvec(o);
corrvec=corrvec(o);

l=length(distvec);

pos=sum(corrvec>=0);
if pos<l
    f=find(corrvec<0,1,'first');
    corrlength=interp1q(corrvec(f:-1:(f-1))',distvec(f:-1:(f-1))',0); 
else corrlength=max(distvec);
end

recrec=sum(sum(corrs(allreceivers,allreceivers)))/2/numreceivers^2;
recnonrec=sum(sum(corrs(allreceivers,~allreceivers)))/numreceivers/(N-numreceivers);
nonrecnonrec=sum(sum(corrs(~allreceivers,~allreceivers)))/2/(N-numreceivers)^2;
total=sum(sum(corrs))/2/N^2;
crosscorrs=[recrec,recnonrec,nonrecnonrec,total];
end

