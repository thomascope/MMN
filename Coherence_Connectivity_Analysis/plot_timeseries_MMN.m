%datapathstem = '/imaging/mlr/users/tc02/Holly_MMN/Coherence_Connectivity_secondfilter';
datapathstem = '/imaging/mlr/users/tc02/Holly_MMN/Coherence_Connectivity_minimumnorm/';
thispath = pwd;
cd(datapathstem)
%load([datapathstem 'groups.mat']);

allfiles = dir(['timeseries*.mat']);

this_con = 0;
this_PNFA = 0;
this_bv = 0;
this_PCA = 0;
All_Ds = [];
All_conDs = [];
All_PNFADs = [];
All_bvDs = [];
All_PCADs = [];
All_Ds_demeaned = [];
All_conDs_demeaned = [];
All_PNFADs_demeaned = [];
All_bvDs_demeaned = [];
All_PCADs_demeaned = [];

for s = 1:size(allfiles,1)
    D = spm_eeg_load(allfiles(s).name);
    D_filt = [];
    
    
%     %Filter below 40Hz
%     fc = 40;
%     Wn = (2/D.fsample)*fc;
%     b = fir1(20,Wn,'low',kaiser(21,3));
%     D_filt = filter(b,1,D(:,:,:));

for this_cond = 1:length(D.condlist)
D_filt(:,:,this_cond) = mean(selectdata(D, [], [], D.condlist{this_cond}),3); %Average by condition
end
%
All_Ds(s,:,:,:) = D_filt(:,:,:);
All_Ds_demeaned(s,:,:,:) = D_filt(:,:,:) / mean(mean(mean(abs(D_filt(:,:,:)),2),3),1);

if ~isempty(strfind(allfiles(1).name,'HC'))
    All_conDs(end+1,:,:,:) = D_filt(:,:,:);
    All_conDs_demeaned(s,:,:,:) = D_filt(:,:,:) / mean(mean(mean(abs(D_filt(:,:,:)),2),3),1);
elseif ~isempty(strfind(allfiles(1).name,'pnfa'))
    All_PNFADs(end+1,:,:,:) = D_filt(:,:,:);
    All_PNFADs_demeaned(s,:,:,:) = D_filt(:,:,:) / mean(mean(mean(abs(D_filt(:,:,:)),2),3),1);
elseif ~isempty(strfind(allfiles(1).name,'bvFTD'))
    All_bvDs(end+1,:,:,:) = D_filt(:,:,:);
    All_bvDs_demeaned(s,:,:,:) = D_filt(:,:,:) / mean(mean(mean(abs(D_filt(:,:,:)),2),3),1);
elseif ~isempty(strfind(allfiles(1).name,'PCA'))
    All_PCADs(end+1,:,:,:) = D_filt(:,:,:);
    All_PCADs_demeaned(s,:,:,:) = D_filt(:,:,:) / mean(mean(mean(abs(D_filt(:,:,:)),2),3),1);
    
end
end
% 
% figure
% plot(D.time,squeeze(mean(mean(All_Ds(:,1,:,[1,3,5]),4),1)),'r--')
% hold on
% plot(D.time,squeeze(mean(mean(All_Ds(:,1,:,[2,4,6]),4),1)),'r:')
% plot(D.time,squeeze(mean(mean(All_Ds(:,2,:,[1,3,5]),4),1)),'b--')
% plot(D.time,squeeze(mean(mean(All_Ds(:,2,:,[2,4,6]),4),1)),'b:')
% 
% figure
% plot(D.time,squeeze(mean(mean(abs(All_Ds(:,1,:,[1,3,5])),4),1)),'r--')
% hold on
% plot(D.time,squeeze(mean(mean(abs(All_Ds(:,1,:,[2,4,6])),4),1)),'r:')
% plot(D.time,squeeze(mean(mean(abs(All_Ds(:,2,:,[1,3,5])),4),1)),'b--')
% plot(D.time,squeeze(mean(mean(abs(All_Ds(:,2,:,[2,4,6])),4),1)),'b:')
% 
% figure
% plot(D.time,squeeze(mean(mean(abs(All_conDs(:,1,:,[1,3,5])),4),1)),'r--')
% hold on
% plot(D.time,squeeze(mean(mean(abs(All_conDs(:,1,:,[2,4,6])),4),1)),'r:')
% plot(D.time,squeeze(mean(mean(abs(All_conDs(:,2,:,[1,3,5])),4),1)),'b--')
% plot(D.time,squeeze(mean(mean(abs(All_conDs(:,2,:,[2,4,6])),4),1)),'b:')
% 
% figure
% plot(D.time,squeeze(mean(mean(abs(All_patDs(:,1,:,[1,3,5])),4),1)),'r--')
% hold on
% plot(D.time,squeeze(mean(mean(abs(All_patDs(:,1,:,[2,4,6])),4),1)),'r:')
% plot(D.time,squeeze(mean(mean(abs(All_patDs(:,2,:,[1,3,5])),4),1)),'b--')
% plot(D.time,squeeze(mean(mean(abs(All_patDs(:,2,:,[2,4,6])),4),1)),'b:')
% 

figure
plot(D.time,squeeze(mean(mean(abs(All_Ds(:,:,:,1)),4),1)))
legend(D.chanlabels)
title('Overall Standard power by source');

figure
plot(D.time,squeeze(mean(mean(abs(All_Ds_demeaned(:,:,:,1)),4),1)))
legend(D.chanlabels)
title('Overall Demeaned Standard power by source');

figure
plot(D.time,squeeze(mean(mean(abs(All_Ds(:,:,:,2)),4),1)))
legend(D.chanlabels)
title('Overall Deviant power by source');

figure
plot(D.time,squeeze(mean(mean(abs(All_Ds_demeaned(:,:,:,2)),4),1)))
legend(D.chanlabels)
title('Overall Demeaned Deviant power by source');

figure
plot(D.time,squeeze(mean(mean(abs(All_Ds(:,:,:,2)),4),1))-squeeze(mean(mean(abs(All_Ds(:,:,:,1)),4),1)))
legend(D.chanlabels)
title('Overall MMN power by source');

figure
plot(D.time,squeeze(mean(mean(abs(All_Ds_demeaned(:,:,:,2)),4),1))-squeeze(mean(mean(abs(All_Ds_demeaned(:,:,:,1)),4),1)))
legend(D.chanlabels)
title('Overall Demeaned MMN power by source');

for i = 1:7
    
figure
plot(D.time,squeeze(mean(mean(abs(All_Ds(:,:,:,i)),4),1)))
legend(D.chanlabels)
title(['Overall ' D.condlist{i} ' power by source']);
end

figure
plot(D.time,squeeze(mean(mean(abs(All_Ds(:,:,:,3:7)),4),1)))
legend(D.chanlabels)
title(['Overall Recalculated Deviant power by source']);


figure
plot(D.time,squeeze(mean(mean(abs(All_Ds_demeaned(:,1,:,1)),4),1)),'r--')
hold on
plot(D.time,squeeze(mean(mean(abs(All_Ds_demeaned(:,1,:,[2,4,6])),4),1)),'r')
plot(D.time,squeeze(mean(mean(abs(All_Ds_demeaned(:,2,:,[1,3,5])),4),1)),'b--')
plot(D.time,squeeze(mean(mean(abs(All_Ds_demeaned(:,2,:,[2,4,6])),4),1)),'b')
ylim([0 1.6]);
title('Overall power by source and condition');
legend({'MisMatch Frontal';'Match Frontal';'MisMatch Temporal';'Match Temporal'});






figure
plot(D.time,squeeze(mean(mean(abs(All_conDs_demeaned(:,1,:,[1,3,5])),4),1)),'r--')
hold on
plot(D.time,squeeze(mean(mean(abs(All_conDs_demeaned(:,1,:,[2,4,6])),4),1)),'r')
plot(D.time,squeeze(mean(mean(abs(All_conDs_demeaned(:,2,:,[1,3,5])),4),1)),'b--')
plot(D.time,squeeze(mean(mean(abs(All_conDs_demeaned(:,2,:,[2,4,6])),4),1)),'b')
ylim([0 1.6]);
title('Control power by source and condition');
legend({'MisMatch Frontal';'Match Frontal';'MisMatch Temporal';'Match Temporal'});

figure
plot(D.time,squeeze(mean(mean(abs(All_patDs_demeaned(:,1,:,[1,3,5])),4),1)),'r--')
hold on
plot(D.time,squeeze(mean(mean(abs(All_patDs_demeaned(:,1,:,[2,4,6])),4),1)),'r')
plot(D.time,squeeze(mean(mean(abs(All_patDs_demeaned(:,2,:,[1,3,5])),4),1)),'b--')
plot(D.time,squeeze(mean(mean(abs(All_patDs_demeaned(:,2,:,[2,4,6])),4),1)),'b')
ylim([0 1.6]);
title('Patient power by source and condition');
legend({'MisMatch Frontal';'Match Frontal';'MisMatch Temporal';'Match Temporal'});

figure
plot(D.time,squeeze(mean(mean(abs(All_conDs_demeaned(:,1,:,:)),4),1)),'r--')
hold on
plot(D.time,squeeze(mean(mean(abs(All_patDs_demeaned(:,1,:,:)),4),1)),'b--')
plot(D.time,squeeze(mean(mean(abs(All_conDs_demeaned(:,2,:,:)),4),1)),'r')
plot(D.time,squeeze(mean(mean(abs(All_patDs_demeaned(:,2,:,:)),4),1)),'b')
title('By group overall power by source');
legend({'Controls Frontal';'Patients Frontal';'Controls Temporal';'Patients Temporal'});

figure
plot(D.time,squeeze(mean(mean(abs(All_conDs_demeaned(:,1,:,[1,3,5])),4),1))-squeeze(mean(mean(abs(All_conDs_demeaned(:,1,:,[2,4,6])),4),1)),'r--')
hold on
plot(D.time,squeeze(mean(mean(abs(All_patDs_demeaned(:,1,:,[1,3,5])),4),1))-squeeze(mean(mean(abs(All_patDs_demeaned(:,1,:,[2,4,6])),4),1)),'b--')
plot(D.time,squeeze(mean(mean(abs(All_conDs_demeaned(:,2,:,[1,3,5])),4),1))-squeeze(mean(mean(abs(All_conDs_demeaned(:,2,:,[2,4,6])),4),1)),'r')
plot(D.time,squeeze(mean(mean(abs(All_patDs_demeaned(:,2,:,[1,3,5])),4),1))-squeeze(mean(mean(abs(All_patDs_demeaned(:,2,:,[2,4,6])),4),1)),'b')
plot(D.time,zeros(size(D.time)),'k--')
title('By group congruency contrast by source');
legend({'Controls Frontal';'Patients Frontal';'Controls Temporal';'Patients Temporal'});

cd(thispath)
