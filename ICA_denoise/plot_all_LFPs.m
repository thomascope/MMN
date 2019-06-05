function plot_all_LFPs(Participant,pathstem,p,prefix)
% A function for plotting my extracted LFPs

Sname = {'left A1';
    'left STG';
    'left IFG';
    'left IPC';
    'right A1';
    'right STG';
    'right IFG';
    'right IPC'};

for ss = 1:length(Participant)
    megpath{ss} = [pathstem Participant{ss}.groupfolder '/' Participant{ss}.name '/' 's_' p.time_wind_path{p.wind_cnt} '_' p.inv_meth{p.inv_cnt} '_' prefix Participant{ss}.name '.mat'];
    diagnosis{ss} = Participant{ss}.diag;
    
    [f1,f2,f3] = fileparts(megpath{ss});
    fn{ss} = sprintf('%s/%s/%dLFP_%s%s',[pathstem 'LFPs'],diagnosis{ss}, length(Sname), f2, f3);
    
end

[groups,~, group_inds] = unique(diagnosis,'stable');

for ss = 1:length(Participant)
    D{ss} = spm_eeg_load(megpath{ss});
    for i = 1:8
       all_STD(i,ss,:)=D{ss}(i,:,1);
       all_DEV(i,ss,:)=D{ss}(i,:,2);
       all_MMN(i,ss,:)=D{ss}(i,:,1)-D{ss}(i,:,2);
       all_abs_MMN(i,ss,:)=abs(D{ss}(i,:,1))-abs(D{ss}(i,:,2));
    end
    LA1_STD(ss,:)=D{ss}(1,:,1);
    LSTG_STD(ss,:)=D{ss}(2,:,1);
    LIFG_STD(ss,:)=D{ss}(3,:,1);
    LIPC_STD(ss,:)=D{ss}(4,:,1);
    RA1_STD(ss,:)=D{ss}(5,:,1);
    RSTG_STD(ss,:)=D{ss}(6,:,1);
    RIFG_STD(ss,:)=D{ss}(7,:,1);
    RIPC_STD(ss,:)=D{ss}(8,:,1);
    LA1_DEV(ss,:)=D{ss}(1,:,2);
    LSTG_DEV(ss,:)=D{ss}(2,:,2);
    LIFG_DEV(ss,:)=D{ss}(3,:,2);
    LIPC_DEV(ss,:)=D{ss}(4,:,2);
    RA1_DEV(ss,:)=D{ss}(5,:,2);
    RSTG_DEV(ss,:)=D{ss}(6,:,2);
    RIFG_DEV(ss,:)=D{ss}(7,:,2);
    RIPC_DEV(ss,:)=D{ss}(8,:,2);
end

cmap = colormap(parula(length(groups)));

% Overall_plot = figure;
% subplot(2,4,1)
% plot(D{1}.time,rms(LA1_STD))
% 
% Brokendown_plot = figure;
% for i = 1:8
% subplot(2,4,i)
% hold on
% linehandle(1) = plot(D{1}.time,rms(squeeze(all_STD(i,1:num_HCs,:))),'k');
% linehandle(2) = plot(D{1}.time,rms(squeeze(all_STD(i,starting_bvFTD:starting_bvFTD+num_bvFTD-1,:))),'r');
% linehandle(3) = plot(D{1}.time,rms(squeeze(all_STD(i,starting_pca:starting_pca+num_pca-1,:))),'b');
% linehandle(4) = plot(D{1}.time,rms(squeeze(all_STD(i,starting_nfvPPA:starting_nfvPPA+num_nfvPPA-1,:))),'g');
% linehandle(5) = plot(D{1}.time,rms(squeeze(all_DEV(i,1:num_HCs,:))),'k--');
% linehandle(6) = plot(D{1}.time,rms(squeeze(all_DEV(i,starting_bvFTD:starting_bvFTD+num_bvFTD-1,:))),'r--');
% linehandle(7) = plot(D{1}.time,rms(squeeze(all_DEV(i,starting_pca:starting_pca+num_pca-1,:))),'b--');
% linehandle(8) = plot(D{1}.time,rms(squeeze(all_DEV(i,starting_nfvPPA:starting_nfvPPA+num_nfvPPA-1,:))),'g--');
% if i==4
%     legend(linehandle,{'HC STD','bvFTD STD','PCA STD','nfvPPA STD','HC DEV','bvFTD DEV','PCA DEV','nfvPPA DEV'})
% end
% end

addpath('./stdshade')
MMN_plot = figure;
for i = 1:8
subplot(2,4,i)
hold on
for grp = 1:length(groups)
linehandle(grp) = stdshade_TEC(squeeze(all_MMN(i,group_inds==grp,:)),0.2,cmap(grp,:),D{1}.time,1,1);
end
title(Sname{i},'FontSize',34)
xlabel('Time (s)')
xlim([-0.1 0.500])
plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
ylabel('Mismatch Response (AU)')
set(gca,'FontWeight','bold')
set(gca,'LineWidth',2)
if i==4
    legend(linehandle,groups)
end
end

STD_plot = figure;
for i = 1:8
subplot(2,4,i)
hold on
for grp = 1:length(groups)
linehandle(grp) = stdshade_TEC(squeeze(all_STD(i,group_inds==grp,:)),0.2,cmap(grp,:),D{1}.time,1,1);
end
title(Sname{i},'FontSize',34)
xlabel('Time (s)')
xlim([-0.1 0.500])
plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
ylabel('Standard Response (AU)')
set(gca,'FontWeight','bold')
set(gca,'LineWidth',2)
if i==4
    legend(linehandle,groups)
end
end


DEV_plot = figure;
for i = 1:8
subplot(2,4,i)
hold on
for grp = 1:length(groups)
linehandle(grp) = stdshade_TEC(squeeze(all_DEV(i,group_inds==grp,:)),0.2,cmap(grp,:),D{1}.time,1,1);
end
title(Sname{i},'FontSize',34)
xlabel('Time (s)')
xlim([-0.1 0.500])
plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
ylabel('Deviant Response (AU)')
set(gca,'FontWeight','bold')
set(gca,'LineWidth',2)
if i==4
    legend(linehandle,groups)
end
end

STD_plot = figure;
for i = 1:8
subplot(2,4,i)
hold on
for grp = 1:length(groups)
linehandle(grp) = stdshade_TEC(squeeze(abs(all_STD(i,group_inds==grp,:))),0.2,cmap(grp,:),D{1}.time,1,1);
end
title(Sname{i},'FontSize',34)
xlabel('Time (s)')
xlim([-0.1 0.500])
plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
ylabel('Abs Standard Response (AU)')
set(gca,'FontWeight','bold')
set(gca,'LineWidth',2)
if i==4
    legend(linehandle,groups)
end
end


DEV_plot = figure;
for i = 1:8
subplot(2,4,i)
hold on
for grp = 1:length(groups)
linehandle(grp) = stdshade_TEC(squeeze(abs(all_DEV(i,group_inds==grp,:))),0.2,cmap(grp,:),D{1}.time,1,1);
end
title(Sname{i},'FontSize',34)
xlabel('Time (s)')
xlim([-0.1 0.500])
plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
ylabel('Abs Deviant Response (AU)')
set(gca,'FontWeight','bold')
set(gca,'LineWidth',2)
if i==4
    legend(linehandle,groups)
end
end


%set(gcf,'color','w')



addpath('./stdshade')
ABS_MMN_plot = figure;
for i = 1:8
subplot(2,4,i)
hold on
for grp = 1:length(groups)
linehandle(grp) = stdshade_TEC(squeeze(all_abs_MMN(i,group_inds==grp,:)),0.2,cmap(grp,:),D{1}.time,1,1);
end
title(Sname{i},'FontSize',34)
xlabel('Time (s)')
xlim([-0.1 0.500])
plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
ylabel('Absolute Mismatch Response (AU)')
set(gca,'FontWeight','bold')
set(gca,'LineWidth',2)
if i==4
    legend(linehandle,groups)
end
end


%pause
