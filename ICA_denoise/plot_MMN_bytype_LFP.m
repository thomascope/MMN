function plot_MMN_bytype_LFP(Participant,pathstem,p,prefix,baselined)
% A function for plotting my extracted LFPs
addpath(['/group/language/data/thomascope/MMN/ICA_denoise/stdshade']);

Sname = {'left A1';
    'left STG';
    'left IFG';
    'left IPC';
    'right A1';
    'right STG';
    'right IFG';
    'right IPC'};

conditions = {'STD','DVT','Loc','Int','Dur','Gap','Freq','Loc_L','Freq_hi','Int_hi','Loc_R','Freq_lo','Int_lo'};

mkdir('./outputfigures/source/')

if ~exist('baselined', 'var')
    baselined = 0;
end

for ss = 1:length(Participant)
    
    megpath{ss} = [pathstem Participant{ss}.groupfolder '/' Participant{ss}.namepostmerge '/' 's_' p.time_wind_path{p.wind_cnt} '_' p.inv_meth{p.inv_cnt} '_' prefix Participant{ss}.namepostmerge '.mat'];
    if ~exist(megpath{ss},'file')
       megpath{ss} = [pathstem Participant{ss}.groupfolder '/' Participant{ss}.namepostmerge '/' 's_' p.time_wind_path{p.wind_cnt} '_' p.inv_meth{p.inv_cnt} '_' prefix Participant{ss}.name{1} '.mat']; 
    end
    diagnosis{ss} = Participant{ss}.diag;
    
    [f1,f2,f3] = fileparts(megpath{ss});
    if baselined == 1
        fn{ss} = sprintf('%s/%s/b%dLFP_%s%s',[pathstem 'LFPs'],diagnosis{ss}, length(Sname), f2, f3);
    else
        fn{ss} = sprintf('%s/%s/%dLFP_%s%s',[pathstem 'LFPs'],diagnosis{ss}, length(Sname), f2, f3);
    end
    
end

[groups,~, group_inds] = unique(diagnosis,'stable');

for ss = 1:length(Participant)
    D{ss} = spm_eeg_load(fn{ss});
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
p_thresh = 0.05;

for j = 1:length(conditions)
    addpath('./stdshade')
    
    DEV_plot = figure(20*j);
    set(gcf,'Position',[100 100 1600 800]);
    set(gcf, 'PaperPositionMode', 'auto');
    for i = 1:8
        subplot(2,4,i)
        hold on
        for grp = 1:length(groups)
            linehandle(grp) = stdshade_TEC(squeeze(all_DEV(i,group_inds==grp,:,j)),0.2,cmap(grp,:),D{1}.time,1,1);
        end
        % Plot statistically significant difference from control (FDR corrected for multiple comparisons across time)
        these_ylims = ylim;
        this_ylim_range = these_ylims(2)-these_ylims(1);
        for grp = 2:length(groups)
            for t = 1:size(D{1}.time,2)
                [h_bytime(t),p_bytime(t)]=ttest2(squeeze(all_DEV(i,group_inds==1,t,j)),squeeze(all_DEV(i,group_inds==grp,t,j)));
            end
            all_fdr_p = mafdr(p_bytime(D{1}.time>0));
            all_positive_times = D{1}.time(D{1}.time>0);
            if any(all_positive_times(all_fdr_p<p_thresh))
                plot(all_positive_times(all_fdr_p<p_thresh),these_ylims(2)+(grp*this_ylim_range/40),'.','Color',cmap(grp,:))
            end
        end
        ylim([these_ylims(1),these_ylims(2)+(this_ylim_range/8)])
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
    saveas(DEV_plot,['./outputfigures/source/' conditions{j} ' Response (AU).png']);
    
    Abs_DEV_plot = figure(400*j);
    set(gcf,'Position',[100 100 1600 800]);
    set(gcf, 'PaperPositionMode', 'auto');
    for i = 1:8
        subplot(2,4,i)
        hold on
        for grp = 1:length(groups)
            linehandle(grp) = stdshade_TEC(squeeze(abs(all_DEV(i,group_inds==grp,:,j))),0.2,cmap(grp,:),D{1}.time,1,1);
        end
        % Plot statistically significant difference from control (FDR corrected for multiple comparisons across time)
        these_ylims = ylim;
        this_ylim_range = these_ylims(2)-these_ylims(1);
        for grp = 2:length(groups)
            for t = 1:size(D{1}.time,2)
                [h_bytime(t),p_bytime(t)]=ttest2(squeeze(abs(all_DEV(i,group_inds==1,t,j))),squeeze(abs(all_DEV(i,group_inds==grp,t,j))));
            end
            all_fdr_p = mafdr(p_bytime(D{1}.time>0));
            all_positive_times = D{1}.time(D{1}.time>0);
            if any(all_positive_times(all_fdr_p<p_thresh))
                plot(all_positive_times(all_fdr_p<p_thresh),these_ylims(2)+(grp*this_ylim_range/40),'.','Color',cmap(grp,:))
            end
        end
        ylim([these_ylims(1),these_ylims(2)+(this_ylim_range/8)])
        title(Sname{i},'FontSize',34)
        xlabel('Time (s)')
        xlim([-0.1 0.500])
        plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
        ylabel([conditions{j} ' Abs Response (AU)'])
        set(gca,'FontWeight','bold')
        set(gca,'LineWidth',2)
        if i==4
            legend(linehandle,groups)
        end
    end
    saveas(Abs_DEV_plot,['./outputfigures/source/' conditions{j} ' Abs Response (AU).png']);
    
    if j > 1 %Don't do MMN for STDs only
        MMN_plot = figure(8000*j);
        set(gcf,'Position',[100 100 1600 800]);
        set(gcf, 'PaperPositionMode', 'auto');
        for i = 1:8
            subplot(2,4,i)
            hold on
            for grp = 1:length(groups)
                linehandle(grp) = stdshade_TEC(squeeze(all_MMN(i,group_inds==grp,:,j)),0.2,cmap(grp,:),D{1}.time,1,1);
            end
            % Plot statistically significant difference from control (FDR corrected for multiple comparisons across time)
            these_ylims = ylim;
            this_ylim_range = these_ylims(2)-these_ylims(1);
            for grp = 2:length(groups)
                for t = 1:size(D{1}.time,2)
                    [h_bytime(t),p_bytime(t)]=ttest2(squeeze(all_MMN(i,group_inds==1,t,j)),squeeze(all_MMN(i,group_inds==grp,t,j)));
                end
                all_fdr_p = mafdr(p_bytime(D{1}.time>0));
                all_positive_times = D{1}.time(D{1}.time>0);
                if any(all_positive_times(all_fdr_p<p_thresh))
                    plot(all_positive_times(all_fdr_p<p_thresh),these_ylims(2)+(grp*this_ylim_range/40),'.','Color',cmap(grp,:))
                end
            end
            ylim([these_ylims(1),these_ylims(2)+(this_ylim_range/8)])
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
        saveas(MMN_plot,['./outputfigures/source/' conditions{j} ' Mismatch Response (AU).png']);
        
        ABS_MMN_plot = figure(160000*j);
        set(gcf,'Position',[100 100 1600 800]);
        set(gcf, 'PaperPositionMode', 'auto');
        for i = 1:8
            subplot(2,4,i)
            hold on
            for grp = 1:length(groups)
                linehandle(grp) = stdshade_TEC(squeeze(all_abs_MMN(i,group_inds==grp,:,j)),0.2,cmap(grp,:),D{1}.time,1,1);
            end
            % Plot statistically significant difference from control (FDR corrected for multiple comparisons across time)
            these_ylims = ylim;
            this_ylim_range = these_ylims(2)-these_ylims(1);
            for grp = 2:length(groups)
                for t = 1:size(D{1}.time,2)
                    [h_bytime(t),p_bytime(t)]=ttest2(squeeze(all_abs_MMN(i,group_inds==1,t,j)),squeeze(all_abs_MMN(i,group_inds==grp,t,j)));
                end
                all_fdr_p = mafdr(p_bytime(D{1}.time>0));
                all_positive_times = D{1}.time(D{1}.time>0);
                if any(all_positive_times(all_fdr_p<p_thresh))
                    plot(all_positive_times(all_fdr_p<p_thresh),these_ylims(2)+(grp*this_ylim_range/40),'.','Color',cmap(grp,:))
                end
            end
            ylim([these_ylims(1),these_ylims(2)+(this_ylim_range/8)])
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
        saveas(ABS_MMN_plot,['./outputfigures/source/' conditions{j} ' Abs Mismatch Response (AU).png']);
    end
    %pause
    close all %To prevent Java memory error
end
%Now plot each region by condition
for i = 1:8
    addpath('./stdshade')
    
    DEV_plot = figure(20*i);
    set(gcf,'Position',[100 100 1600 800]);
    set(gcf, 'PaperPositionMode', 'auto');
    for j = 1:7
        subplot(2,4,j)
        hold on
        for grp = 1:length(groups)
            linehandle(grp) = stdshade_TEC(squeeze(all_DEV(i,group_inds==grp,:,j)),0.2,cmap(grp,:),D{1}.time,1,1);
        end
        % Plot statistically significant difference from control (FDR corrected for multiple comparisons across time)
        these_ylims = ylim;
        this_ylim_range = these_ylims(2)-these_ylims(1);
        for grp = 2:length(groups)
            for t = 1:size(D{1}.time,2)
                [h_bytime(t),p_bytime(t)]=ttest2(squeeze(all_DEV(i,group_inds==1,t,j)),squeeze(all_DEV(i,group_inds==grp,t,j)));
            end
            all_fdr_p = mafdr(p_bytime(D{1}.time>0));
            all_positive_times = D{1}.time(D{1}.time>0);
            if any(all_positive_times(all_fdr_p<p_thresh))
                plot(all_positive_times(all_fdr_p<p_thresh),these_ylims(2)+(grp*this_ylim_range/40),'.','Color',cmap(grp,:))
            end
        end
        ylim([these_ylims(1),these_ylims(2)+(this_ylim_range/8)])
        title(conditions{j},'FontSize',34)
        xlabel('Time (s)')
        xlim([-0.1 0.500])
        plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
        ylabel([conditions{j} ' Response (AU)'])
        set(gca,'FontWeight','bold')
        set(gca,'LineWidth',2)
    end
    for j = 8
        subplot(2,4,j)
        hold on
        for grp = 1:length(groups)
            linehandle(grp) = stdshade_TEC(squeeze(all_DEV(i,group_inds==grp,:,j)),0,cmap(grp,:),D{1}.time,1,1);
        end
        axis off
        legend(linehandle,groups)
        drawnow
        delete(linehandle)
        drawnow
    end
    saveas(DEV_plot,['./outputfigures/source/' Sname{i} ' Response (AU).png']);
    
    Abs_DEV_plot = figure(400*i);
    set(gcf,'Position',[100 100 1600 800]);
    set(gcf, 'PaperPositionMode', 'auto');
    for j = 1:7
        subplot(2,4,j)
        hold on
        for grp = 1:length(groups)
            linehandle(grp) = stdshade_TEC(squeeze(abs(all_DEV(i,group_inds==grp,:,j))),0.2,cmap(grp,:),D{1}.time,1,1);
        end
        % Plot statistically significant difference from control (FDR corrected for multiple comparisons across time)
        these_ylims = ylim;
        this_ylim_range = these_ylims(2)-these_ylims(1);
        for grp = 2:length(groups)
            for t = 1:size(D{1}.time,2)
                [h_bytime(t),p_bytime(t)]=ttest2(squeeze(abs(all_DEV(i,group_inds==1,t,j))),squeeze(abs(all_DEV(i,group_inds==grp,t,j))));
            end
            all_fdr_p = mafdr(p_bytime(D{1}.time>0));
            all_positive_times = D{1}.time(D{1}.time>0);
            if any(all_positive_times(all_fdr_p<p_thresh))
                plot(all_positive_times(all_fdr_p<p_thresh),these_ylims(2)+(grp*this_ylim_range/40),'.','Color',cmap(grp,:))
            end
        end
        ylim([these_ylims(1),these_ylims(2)+(this_ylim_range/8)])
        title(['Abs ' conditions{j}],'FontSize',28)
        xlabel('Time (s)')
        xlim([-0.1 0.500])
        plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
        ylabel([conditions{j} ' Abs Response (AU)'])
        set(gca,'FontWeight','bold')
        set(gca,'LineWidth',2)
    end
    for j = 8
        subplot(2,4,j)
        hold on
        for grp = 1:length(groups)
            linehandle(grp) = stdshade_TEC(squeeze(all_DEV(i,group_inds==grp,:,j)),0,cmap(grp,:),D{1}.time,1,1);
        end
        axis off
        legend(linehandle,groups)
        drawnow
        delete(linehandle)
        drawnow
    end
    saveas(Abs_DEV_plot,['./outputfigures/source/' Sname{i} ' Abs Response (AU).png']);
    
        MMN_plot = figure(8000*i);
        set(gcf,'Position',[100 100 1600 800]);
        set(gcf, 'PaperPositionMode', 'auto');
        for j = 1:7
            subplot(2,4,j)
            hold on
            for grp = 1:length(groups)
                linehandle(grp) = stdshade_TEC(squeeze(all_MMN(i,group_inds==grp,:,j)),0.2,cmap(grp,:),D{1}.time,1,1);
            end
            % Plot statistically significant difference from control (FDR corrected for multiple comparisons across time)
            these_ylims = ylim;
            this_ylim_range = these_ylims(2)-these_ylims(1);
            if j >= 2
            for grp = 2:length(groups)
                for t = 1:size(D{1}.time,2)
                    [h_bytime(t),p_bytime(t)]=ttest2(squeeze(all_MMN(i,group_inds==1,t,j)),squeeze(all_MMN(i,group_inds==grp,t,j)));
                end
                all_fdr_p = mafdr(p_bytime(D{1}.time>0));
                all_positive_times = D{1}.time(D{1}.time>0);
                if any(all_positive_times(all_fdr_p<p_thresh))
                    plot(all_positive_times(all_fdr_p<p_thresh),these_ylims(2)+(grp*this_ylim_range/40),'.','Color',cmap(grp,:))
                end
            end
            end
            ylim([these_ylims(1),these_ylims(2)+(this_ylim_range/8)])
            title(['MMN ' conditions{j}],'FontSize',28)
            xlabel('Time (s)')
            xlim([-0.1 0.500])
            plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
            ylabel([conditions{j} ' Mismatch Response (AU)'])
            set(gca,'FontWeight','bold')
            set(gca,'LineWidth',2)
        end
        for j = 8
            subplot(2,4,j)
            hold on
            for grp = 1:length(groups)
                linehandle(grp) = stdshade_TEC(squeeze(all_MMN(i,group_inds==grp,:,j)),0,cmap(grp,:),D{1}.time,1,1);
            end
            axis off
            legend(linehandle,groups)
            drawnow
            delete(linehandle)
            drawnow
        end
        saveas(MMN_plot,['./outputfigures/source/' Sname{i} ' Mismatch Response (AU).png']);
        
        ABS_MMN_plot = figure(160000*i);
        set(gcf,'Position',[100 100 1600 800]);
        set(gcf, 'PaperPositionMode', 'auto');
        for j = 1:7
            subplot(2,4,j)
            hold on
            for grp = 1:length(groups)
                linehandle(grp) = stdshade_TEC(squeeze(all_abs_MMN(i,group_inds==grp,:,j)),0.2,cmap(grp,:),D{1}.time,1,1);
            end
            % Plot statistically significant difference from control (FDR corrected for multiple comparisons across time)
            these_ylims = ylim;
            this_ylim_range = these_ylims(2)-these_ylims(1);
            if j >= 2
            for grp = 2:length(groups)
                for t = 1:size(D{1}.time,2)
                    [h_bytime(t),p_bytime(t)]=ttest2(squeeze(all_abs_MMN(i,group_inds==1,t,j)),squeeze(all_abs_MMN(i,group_inds==grp,t,j)));
                end
                all_fdr_p = mafdr(p_bytime(D{1}.time>0));
                all_positive_times = D{1}.time(D{1}.time>0);
                if any(all_positive_times(all_fdr_p<p_thresh))
                    plot(all_positive_times(all_fdr_p<p_thresh),these_ylims(2)+(grp*this_ylim_range/40),'.','Color',cmap(grp,:))
                end
            end
            end
            ylim([these_ylims(1),these_ylims(2)+(this_ylim_range/8)])
            title(['Abs MMN ' conditions{j}],'FontSize',22)
            xlabel('Time (s)')
            xlim([-0.1 0.500])
            plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
            ylabel([conditions{j} ' Abs Mismatch Response (AU)'])
            set(gca,'FontWeight','bold')
            set(gca,'LineWidth',2)
        end
        for j = 8
            subplot(2,4,j)
            hold on
            for grp = 1:length(groups)
                linehandle(grp) = stdshade_TEC(squeeze(all_abs_MMN(i,group_inds==grp,:,j)),0,cmap(grp,:),D{1}.time,1,1);
            end
            axis off
            legend(linehandle,groups)
            drawnow
            delete(linehandle)
            drawnow
        end
        saveas(ABS_MMN_plot,['./outputfigures/source/' Sname{i} ' Abs Mismatch Response (AU).png']);

    %pause
    close all %To prevent Java memory error
end
%pause
