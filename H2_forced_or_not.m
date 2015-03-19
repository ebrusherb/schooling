Q=zeros(N-1,N);
Q(1,1:2)=[1/sqrt(2) -1/sqrt(2)];
for j=2:(N-1)
   v=rand(1,N);
   v=v-sum(v)/N;
   for k=1:(j-1)
      v=v-(v*transpose(Q(k,:)))*Q(k,:); 
   end
   v=v/norm(v,2);
   Q(j,:)=v;
end
Qfull=[Q; 1/sqrt(N)*ones(1,N)];

N=20;
noise=eye(N);
positions=unifrnd(0,1,N,2);
    d=squareform(pdist(positions));

    A=zeros(N);
    Abar=zeros(N);
    for ind=1:N
        [~, order]=sort(d(ind,:));
        neighbors=order(2:strategy(ind)+1);
        A(ind,neighbors)=1;
        Abar(ind,neighbors)=1/strategy(ind);
    end
  
    L=Abar-eye(N);
    Ltilde=Q*L*transpose(Q);
    
     sigma=lyap(Ltilde,Q*noise*transpose(Q));
        H=sqrt(trace(sigma));
        
%         receivers=randsample(N,numsigs_permove,'true');
        

%             receiver=receivers(1);
receiver=1;
            allreceivers=d(receiver,:)<=radius;
            beta=zeros(N,1);
            beta(allreceivers)=b;
            B=diag(beta);
            Lf=L-B;
            
            sigma_forced=lyap(Lf,noise);
                H_forced=trace(sigma_forced);
                
%                 sigma_forced2=lyap(Q*Lf*transpose(Q),Q*noise*transpose(Q));
            sigmax=sigma_forced+inv(Lf)*B*ones(N,1)*ones(1,N)*transpose(B)*transpose(inv(Lf));
            sigmaz=Qfull*sigmax*transpose(Qfull);
            sigmay=sigmaz(1:end-1,1:end-1);
            sigmaz2=lyap(Qfull*Lf*inv(Qfull),+Qfull*transpose(Qfull));
            sigmax2=transpose(Qfull)*sigmaz2*Qfull;
%             sigmay2=
                
                subplot(2,2,1);imagesc(sigma);colorbar;
                subplot(2,2,2);imagesc(sigmay);colorbar;
                subplot(2,2,3);imagesc(sigmax);colorbar;
                subplot(2,2,4);imagesc(sigmax2);colorbar
%                 subplot(2,2,2);imagesc(sigma_forced2);colorbar;
%                 
% %                 subplot(2,2,3);imagesc(transpose(Q)*sigma*Q);colorbar;
%                 subplot(2,2,3);imagesc((Q)*sigma_forced*transpose(Q));colorbar;
% %                 subplot(2,2,3);imagesc(transpose(Q)*sigma_forced2*Q);colorbar;
%                 subplot(2,2,4);imagesc(sigma_forced);colorbar
%                 
                %%
                M=transpose(Q)*sigma_forced2*Q;
               subplot(1,2,1);imagesc(sigma_forced2);colorbar;
                subplot(1,2,2);imagesc(sigma_forced(2:end,2:end));colorbar; 
                
                %%
                M=rand(10,10);
                S=lyap(Q*M*transpose(Q),Q*noise*transpose(Q));
                subplot(1,2,1);imagesc(Q*M*transpose(Q)*S+S*Q*transpose(M)*transpose(Q));
                colorbar;
                subplot(1,2,2);imagesc(-Q*noise*transpose(Q));
                colorbar
                
                %%
                subplot(1,2,1)
                imagesc(transpose(Q)*sigma_forced2)
                colorbar
                subplot(1,2,2)
                imagesc(C*transpose(Q))
                colorbar
                
                %%
Qfull=[Q;1/sqrt(N)*ones(1,N)];
Ltilde=Qfull*Lf*inv(Qfull);
Bbar=Qfull*B*ones(N,1);

sigma=lyap(Ltilde,-inv(Ltilde)*Bbar*transpose(Bbar)-Bbar*transpose(Bbar)*transpose(inv(Ltilde))+Qfull*transpose(Qfull));
C=lyap(Lf,-inv(Lf)*B*ones(N,1)*ones(1,N)*transpose(B)-B*ones(N,1)*ones(1,N)*transpose(B)*transpose(inv(Lf))+eye(N));
C2=lyap(Lf,eye(N));
C3=C-inv(Lf)*B*ones(N,1)*ones(1,N)*transpose(B)*transpose(inv(Lf));
Cbar=Qfull*C*transpose(Qfull);
% subplot(1,2,1)
% imagesc(Lf*C2+C2*transpose(Lf))
% colorbar
% subplot(1,2,2)
% imagesc((eye(N)))
% colorbar
% subplot(1,2,1)
% subplot(1,2,1)
% imagesc(Lf*C+C*transpose(Lf)+inv(Lf)*B*ones(N,1)*ones(1,N)*transpose(B)+B*ones(N,1)*ones(1,N)*transpose(B)*transpose(inv(Lf)))
% colorbar
% subplot(1,2,2)
% imagesc(-(-eye(N)))
% colorbar
% subplot(1,2,1)
% imagesc(Lf*C3+C3*transpose(Lf))
% colorbar
% subplot(1,2,2)
% imagesc(Lf*C+C*transpose(Lf)-inv(Lf)*B*ones(N,1)*ones(1,N)*transpose(B)-B*ones(N,1)*ones(1,N)*transpose(B)*transpose(inv(Lf)))
% colorbar
subplot(1,2,1)
imagesc(C2(1:end,1:end))
colorbar
subplot(1,2,2)
imagesc(C3(1:end,1:end))
colorbar
% subplot(1,2,1)
% imagesc(sigma(1:end-1,1:end-1))
% colorbar
% subplot(1,2,2)
% imagesc(Cbar(1:end-1,1:end-1))
% colorbar
%%
alpha=.1;
Ltilde2=Qfull*L*transpose(Qfull);
Lfbar=Qfull*Lf*transpose(Qfull);
Bbar=alpha*Qfull*B;

sigmay=lyap(Ltilde,Q*transpose(Q));
sigmaz=lyap(Ltilde2,Qfull*transpose(Qfull));
sigmax=transpose(Qfull)*sigmaz*Qfull;
subplot(2,2,1)
imagesc(sigmay);colorbar
subplot(2,2,2)
imagesc(sigmaz(1:end-1,1:end-1));colorbar

sigmaz2=lyap(Lfbar,-inv(Lfbar)*Bbar*ones(N,1)*ones(1,N)*transpose(Bbar)-Bbar*ones(N,1)*ones(1,N)*transpose(Bbar)*transpose(inv(Lfbar))+Qfull*transpose(Qfull));
sigmay2=sigmaz2(1:end-1,1:end-1);
sigmax=lyap(Lf,-inv(Lf)*B*ones(N,1)*ones(1,N)*B-B*ones(N,1)*ones(1,N)*B*transpose(inv(Lf))+eye(N));
sigmaz3=Qfull*sigmax*transpose(Qfull);
sigmax2=transpose(Qfull)*sigmaz2*Qfull;
subplot(2,2,3)
imagesc(sigmax)
colorbar
subplot(2,2,4)
imagesc(sigmax2);
colorbar
% subplot(2,2,3)
% imagesc(sigmax)
% colorbar
% subplot(2,2,4)
% imagesc(sigmay2);
% colorbar
% subplot(2,2,3)
% imagesc(sigmaz2(1:end-1,1:end-1))
% colorbar
% subplot(2,2,4)
% imagesc(sigmaz3(1:end-1,1:end-1));
% colorbar