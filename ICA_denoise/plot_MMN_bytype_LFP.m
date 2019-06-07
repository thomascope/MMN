function plot_MMN_bytype_LFP(Participant,pathstem,p,prefix)
% A function for plotting my extracted LFPs

Sname = {'left A1';
    'left STG';
    'left IFG';
    'left IPC';
    'right A1';
    'right STG';
    'right IFG';
    'right IPC'};

conditions = {'STD','DVT','Loc','Int','Dur','Gap','Freq','Loc_L','Freq_hi','Int_hi','Loc_R','Freq_lo','Int_lo'};

for ss = 1:length(Participant)
    try
        Participant{ss}.name = Participant{ss}.namepostmerge;
    end
    megpath{ss} = [pathstem Participant{ss}.groupfolder '/' Participant{ss}.name '/' 's_' p.time_wind_path{p.wind_cnt} '_' p.inv_meth{p.inv_cnt} '_' prefix Participant{ss}.name '.mat'];
    diagnosis{ss} = Participant{ss}.diag;
    
    [f1,f2,f3] = fileparts(megpath{ss});
    fn{ss} = sprintf('%s/%s/%dLFP_%s%s',[pathstem 'LFPs'],diagnosis{ss}, length(Sname), f2, f3);
    
end

[groups,~, group_inds] = unique(diagnosis,'stable');

for ss = 1:length(Participant)
    D{ss} = spm_eeg_load(megpath{ss});
    for i = 1:length(Sname)
        all_STD(i,ss,:)=D{ss}(i,:,1);
        for j = 1:length(conditions)
            all_DEV(i,ss,:,j)=D{ss}(i,:,j);
            all_MMN(i,ss,:,j)=D{ss}(i,:,1)-D{ss}(i,:,j);
            all_abs_MMN(i,ss,:,j)=abs(D{ss}(i,:,1))-abs(D{ss}(i,:,j));
        end
    end
end

cmap = colormap(parula(length(groups)));

for j = 1:length(conditions)
    addpath('./stdshade')
    
    DEV_plot = figure(20*j);
    for i = 1:8
        subplot(2,4,i)
        hold on
        for grp = 1:length(groups)
            linehandle(grp) = stdshade_TEC(squeeze(all_DEV(i,group_inds==grp,:,j)),0.2,cmap(grp,:),D{1}.time,1,1);
        end
        title(Sname{i},'FontSize',34)
        xlabel('Time (s)')
        xlim([-0.1 0.500])
        plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
        ylabel([conditions{j} ' Response (AU)'])
        set(gca,'FontWeight','bold')
        set(gca,'LineWidth',2)
        if i==4
            legend(linehandle,groups)
        end
    end
    
    MMN_plot = figure(j);
    for i = 1:8
        subplot(2,4,i)
        hold on
        for grp = 1:length(groups)
            linehandle(grp) = stdshade_TEC(squeeze(all_MMN(i,group_inds==grp,:,j)),0.2,cmap(grp,:),D{1}.time,1,1);
        end
        title(Sname{i},'FontSize',34)
        xlabel('Time (s)')
        xlim([-0.1 0.500])
        plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
        ylabel([conditions{j} ' Mismatch Response (AU)'])
        set(gca,'FontWeight','bold')
        set(gca,'LineWidth',2)
        if i==4
            legend(linehandle,groups)
        end
    end
    
    ABS_MMN_plot = figure(400*j);
    for i = 1:8
        subplot(2,4,i)
        hold on
        for grp = 1:length(groups)
            linehandle(grp) = stdshade_TEC(squeeze(all_abs_MMN(i,group_inds==grp,:,j)),0.2,cmap(grp,:),D{1}.time,1,1);
        end
        title(Sname{i},'FontSize',34)
        xlabel('Time (s)')
        xlim([-0.1 0.500])
        plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
        ylabel([conditions{j} ' Abs Mismatch Response (AU)'])
        set(gca,'FontWeight','bold')
        set(gca,'LineWidth',2)
        if i==4
            legend(linehandle,groups)
        end
    end
    pause
    close all %To prevent Java memory error
end
%pause
