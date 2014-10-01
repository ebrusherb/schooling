function t =parallel_example

t0=tic;
parfor idx=1:16
    A(idx)=idx;
    pause(2)
end
t=toc(t0);

% 
% its=20;
% 
% N=40;
% 
% stratvals=[10 15 20];
% bvals=[0 0.01 .1:.1:1 1.5];
% radvals=[.01 .1 .5 1];
% 
% Ns=length(stratvals);
% Nb=length(bvals);
% Nr=length(radvals);
% 
% corrmat=zeros(its,Ns,Nr);
% corrcorrmat=zeros(its,Ns,Nr,Nb);
% H2mat=zeros(its,Ns);
% 
% for itnow=1:its
%     positions=unifrnd(0,1,N,2);
%     d=squareform(pdist(positions));
% 
%     for i=1:Ns
%         strat=stratvals(i);
%         strategy=strat*ones(N,1);
% 
%         M=zeros(N);
%         for ind=1:N
%             [~, order]=sort(d(ind,:));
%             neighbors=order(2:strategy(ind)+1);
%             M(ind,neighbors)=1/strategy(ind);
%         end
%         M(1:N+1:end)=-1;
%         H2mat(itnow,i)=H2norm(M);
%         
%         for k=1:Nr
%             radius=radvals(k);
%             corrcorrvec=zeros(Nb,1);
%             [~,endavgcorr,corrlength] = correlationlength_mat(M,d,1000,radius);
%             endstop=find(isnan(endavgcorr),1,'first')-1;
%             corrmat(itnow,i,k)=corrlength;
%             parfor j=1:Nb 
%                 b=bvals(j);
%                 [~,avgcorr,corrlength] = correlationlength_mat(M,d,b,radius);
%                 nowstop=find(isnan(avgcorr),1,'first')-1;
%                 stop=min(endstop,nowstop);
%                 corrcorrvec(j)=corr(avgcorr(1:stop)',endavgcorr(1:stop)'); %#ok<PFBNS>
%             end
%             corrcorrmat(itnow,i,k,:)=corrcorrvec;
%         end
%     end
% end 
