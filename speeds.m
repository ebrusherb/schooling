function rates=speeds(M,lambdaorv,meanormax)
if nargin==1
    lambdaorv='both';
    meanormax='mean';
end
s=size(M);
N=s(2);
L=lap(M);
delta=zeros(N);
if var(sum(transpose(M)))==0 %if everyone's the same don't bother looping over all receivers
        [path_dists,~,~]=graphshortestpath(sparse(transpose(M)),1);
        not_children=(path_dists==Inf);
        children=(path_dists~=Inf);
        children(1)=0;
        delta(not_children,1)=0;
        if sum(children>0)
            W=L(children,children);
            [v d]=eig(W);
            [lambda,index] = min(real(diag(d)));
            v=v(:,index);
            v=v/mean(v);
            v=abs(v);
            if strcmp(lambdaorv,'both')==1
            %delta(children,i)=exp(-minValue)*abs(v(:,index));
            delta(children,1)=v/sqrt(2*lambda);
            elseif strcmp(lambdaorv,'lambda')==1
            %delta(children,i)=exp(-minValue);
            delta(children,1)=1/sqrt(2*lambda);
            else delta(children,1)=v;
            end
        end
    if strcmp(meanormax,'mean')==1
        rates=mean(delta(:,1))*ones(N,1);
    else 
        rates=max(delta(:,1))*ones(N,1);
    end
end 
if var(sum(transpose(M)))~=0
    for i=1:N
        [path_dists,~,~]=graphshortestpath(sparse(transpose(M)),i);
        not_children=(path_dists==Inf);
        children=(path_dists~=Inf);
        children(i)=0;
        delta(not_children,i)=0;
        if sum(children>0)
            W=L(children,children);
            [v d]=eig(W);
            [lambda,index] = min(real(diag(d)));
            v=v(:,index);
            v=v/mean(v);
            v=abs(v);
            if strcmp(lambdaorv,'both')==1
            %delta(children,i)=exp(-minValue)*abs(v(:,index));
            delta(children,i)=v/sqrt(2*lambda);
            elseif strcmp(lambdaorv,'lambda')==1
            %delta(children,i)=exp(-minValue);
            delta(children,i)=1/sqrt(2*lambda);
            else delta(children,i)=v;
            end
        end
    end
    if strcmp(meanormax,'mean')==1
        rates=mean(delta,2);
    else 
        rates=max(transpose(delta));
        rates=transpose(rates);
    end
end   
