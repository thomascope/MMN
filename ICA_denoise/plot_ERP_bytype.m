function plot_ERP_bytype(Participant,pathstem,p,prefix)
% A function for plotting my extracted LFPs
addpath(['/group/language/data/thomascope/MMN/ICA_denoise/stdshade']);

conditions = {'STD','DVT','Loc','Int','Dur','Gap','Freq','Loc_L','Freq_hi','Int_hi','Loc_R','Freq_lo','Int_lo'};

for ss = 1:length(Participant)
    try
        Participant{ss}.name = Participant{ss}.namepostmerge;
    end
    megpath{ss} = [pathstem Participant{ss}.groupfolder '/' Participant{ss}.name '/' prefix Participant{ss}.name '.mat'];
    diagnosis{ss} = Participant{ss}.diag;
   
end

[groups,~, group_inds] = unique(diagnosis,'stable');

for ss = 1:length(Participant)
    D{ss} = spm_eeg_load(megpath{ss});
    all_STD(ss,:)=mean(D{ss}(D{ss}.selectchannels('MEGCOMB'),:,1),1);
    for j = 1:length(conditions)
        all_DEV(ss,:,j)=mean(D{ss}(D{ss}.selectchannels('MEGCOMB'),:,j),1);
        all_MMN(ss,:,j)=mean(D{ss}(D{ss}.selectchannels('MEGCOMB'),:,1),1)-mean(D{ss}(D{ss}.selectchannels('MEGCOMB'),:,j),1);
        all_abs_MMN(ss,:,j)=mean(abs(D{ss}(D{ss}.selectchannels('MEGCOMB'),:,1)),1)-mean(abs(D{ss}(D{ss}.selectchannels('MEGCOMB'),:,j)),1);
    end
end

cmap = colormap(parula(length(groups)));

for j = 1:length(conditions)
    addpath('./stdshade')
    
    DEV_plot = figure(20*j);

        hold on
        for grp = 1:length(groups)
            linehandle(grp) = stdshade_TEC(squeeze(all_DEV(group_inds==grp,:,j)),0.2,cmap(grp,:),D{1}.time,1,1);
        end
        xlabel('Time (s)')
        xlim([-0.1 0.500])
        plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
        ylabel([conditions{j} ' Response (AU)'])
        set(gca,'FontWeight','bold')
        set(gca,'LineWidth',2)
            legend(linehandle,groups)
    
    MMN_plot = figure(j);

        hold on
        for grp = 1:length(groups)
            linehandle(grp) = stdshade_TEC(squeeze(all_MMN(group_inds==grp,:,j)),0.2,cmap(grp,:),D{1}.time,1,1);
        end
        xlabel('Time (s)')
        xlim([-0.1 0.500])
        plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
        ylabel([conditions{j} ' Mismatch Response (AU)'])
        set(gca,'FontWeight','bold')
        set(gca,'LineWidth',2)
            legend(linehandle,groups)
    
    ABS_MMN_plot = figure(400*j);
        hold on
        for grp = 1:length(groups)
            linehandle(grp) = stdshade_TEC(squeeze(all_abs_MMN(group_inds==grp,:,j)),0.2,cmap(grp,:),D{1}.time,1,1);
        end
        xlabel('Time (s)')
        xlim([-0.1 0.500])
        plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
        ylabel([conditions{j} ' Abs Mismatch Response (AU)'])
        set(gca,'FontWeight','bold')
        set(gca,'LineWidth',2)
            legend(linehandle,groups)
            
            DEV_indiv_plot = figure(8000*j);

        hold on
        for grp = 1:length(groups)
            linehandles{grp} = plot(D{1}.time,squeeze(all_DEV(group_inds==grp,:,j)),'Color',cmap(grp,:));
        end
        xlabel('Time (s)')
        xlim([-0.1 0.500])
        ylabel([conditions{j} ' Individual Response (AU)'])
        set(gca,'FontWeight','bold')
        set(gca,'LineWidth',2)
        for i = 1:length(linehandles)
            linehandle(i) = linehandles{i}(1);
        end
            legend(linehandle,groups)
    pause
    close all %To prevent Java memory error
end
%pause
