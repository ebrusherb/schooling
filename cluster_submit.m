configCluster
%%
dellacluster=parcluster;

ClusterInfo.setEmailAddress('brush@princeton.edu');

localcluster=parcluster('local');

%%
ClusterInfo.setWallTime('00:01:00');
testjob=dellacluster.batch('parallel_example','matlabpool',1);

%%
ClusterInfo.setWallTime('30:00:00');
ESSjob = dellacluster.batch('ESS_parallel','matlabpool',30);

%%
ESSoutputs=ESSjob.fetchOutputs{1};
strats=ESSoutputs.strats;
fitnesseaten=ESSoutputs.fitnesseaten;
rhoeaten=ESSoutputs.rhoeaten;
ESSseaten=ESSoutputs.ESSseaten;
fitnessgettoeat=ESSoutputs.fitnessgettoeat;
rhogettoeat=ESSoutputs.rhogettoeat;
ESSsgettoeat=ESSoutputs.ESSsgettoeat;
%%
ClusterInfo.setWallTime('15:00:00');
annealingjob = dellacluster.batch('simulated_annealing','matlabpool',30);

%%
annealingoutputs=annealingjob.fetchOutputs{1};

%%
% ClusterInfo.setWallTime('06:00:00');
corrlengthjob=dellacluster.batch('corrlength_parallel','matlabpool',16);

%%
corrlength_outputs=corrlengthjob.fetchOutputs{1};
corrlengthmat=corrlength_outputs.corrlengthmat;
similaritymat=corrlength_outputs.similaritymat;
minbmat=corrlength_outputs.minbmat;
H2mat=corrlength_outputs.H2mat;

% Nb2=corrlength_outputs.Nb2;
% storeddistbins=corrlength_outputs.storeddistbins;
% storedavgcorr=corrlength_outputs.storedavgcorr;
% ref=corrlength_outputs.ref;
% c=colormap(jet(Nb2));
% figure
% for i=1:2
%     subplot(1,2,i)
%     hold on
%     for j=1:Nb2
%         plot(storeddistbins{j,i},storedavgcorr{j,i},'Color',c(j,:),'LineWidth',2)
%     end
%     plot(storeddistbins{Nb2,i},ref{i},'k','LineWidth',2)
%     plot([0 1.5],zeros(2,1))
% end
% 
% compcorrmat=zeros(2,Nb2,2);
% for i=1:2
%     for j=1:Nb2
%         compcorrmat(i,j,1)=corr(storedavgcorr{1,i}(1:50)',storedavgcorr{j,i}(1:50)');
%         compcorrmat(i,j,2)=corr(ref{i}(1:50)',storedavgcorr{j,i}(1:50)');
%     end
% end