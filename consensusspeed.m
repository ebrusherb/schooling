function s=consensusspeed(M)
        L=lap(M);
        [~,d]=eig(L);
        d=diag(d);
        good=1-(abs(d)<0.0001);
        d=d(logical(good));
        lambda=min(real(d));
        s=1/sqrt(2*lambda);