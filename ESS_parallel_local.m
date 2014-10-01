t0=tic;
numsigs_permove=1;
nummoves=100;
N=20;
b=1;

strats=2:2:(N-1);
L=length(strats);

Tvals=1;
Nt=length(Tvals);

% radvals=[.0 .1 .5 1.5];
radvals=.1;
Nr=length(radvals);

ESSseaten=cell(Nt,Nr);
ESSsgettoeat=cell(Nt,Nr);

for q=1:Nt
    for p=1:Nr
        T=Tvals(q);
        radius=radvals(p);

        fitnesseaten=zeros(L,L);
        rhoeaten=zeros(L,L);

        fitnessgettoeat=zeros(L,L);
        rhogettoeat=zeros(L,L);

        featen=zeros(L,L,N-1);
        geaten=zeros(L,L,N-1);

        fgettoeat=zeros(L,L,N-1);
        ggettoeat=zeros(L,L,N-1);
        parfor ind=1:(L*L*(N-1))
            [u,v,k]=ind2sub([L,L,N-1],ind);
            squareind=sub2ind([L,L],u,v);
            resident=strats(u);
            invader=strats(v);

            strategy=resident*ones(1,N);
            strategy(1:k)=invader;
            [meanscore, scorevar, probeaten, probgettoeat]=signalingevents(strategy,numsigs_permove,nummoves,radius,b,T);

            perfeaten=1-probeaten;
            featen(ind)=mean(perfeaten(1:k));
            geaten(ind)=mean(perfeaten((k+1):N));

            perfgettoeat=probgettoeat;
            fgettoeat(ind)=mean(perfgettoeat(1:k));
            ggettoeat(ind)=mean(perfgettoeat((k+1):N));

        end


        for ind=1:L*L
            [u,v]=ind2sub([L,L],ind);
            where=zeros(N-1,1);
            for k=1:N-1
                where(k)=sub2ind([L,L,N-1],u,v,k);
            end
            featen_now=featen(where);
            geaten_now=geaten(where);

            fgettoeat_now=fgettoeat(where);
            ggettoeat_now=ggettoeat(where);

            fitnesseaten(ind)=featen_now(1)/geaten_now(1);
            fitnessgettoeat(ind)=fgettoeat_now(1)/ggettoeat_now(1);

            ratio=geaten_now./featen_now;
            tosum=ones(1,N-1);
            for k=1:(N-1)
                tosum(k)=prod(ratio(1:k));
            end
            rhoeaten(ind)=1/(1+sum(tosum));
 
            ratio=ggettoeat_now./fgettoeat_now;
            tosum=ones(1,N-1);
            for k=1:(N-1)
                tosum(k)=prod(ratio(1:k));
            end
            rhogettoeat(ind)=1/(1+sum(tosum));
        end   

        eqstratseaten=eq_strats(N,fitnesseaten,rhoeaten);
        ESSseaten{q,p}=eqstratseaten{2};

        eqstratsgettoeat=eq_strats(N,fitnessgettoeat,rhogettoeat);
        ESSsgettoeat{q,p}=eqstratsgettoeat{2};
    end
end
t=toc(t0);
save('/Users/eleanorbrush/Desktop/ESSoutput.mat','strats','radvals','fitnesseaten','rhoeaten','ESSseaten','fitnessgettoeat','rhogettoeat','ESSsgettoeat','t','t0')
% save('/home/brush/MdcsDataLocation/princeton/della/R2013b/remote/ESSoutput.mat','strats','radvals','fitnesseaten','rhoeaten','ESSseaten','fitnessgettoeat','rhogettoeat','ESSsgettoeat');

