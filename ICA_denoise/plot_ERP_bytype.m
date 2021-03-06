function plot_ERP_bytype(Participant,pathstem,p,prefix,varargin)
% A function for plotting my extracted LFPs
addpath(['/group/language/data/thomascope/MMN/ICA_denoise/stdshade']);

conditions = {'STD','DVT','Loc','Int','Dur','Gap','Freq','Loc_L','Freq_hi','Int_hi','Loc_R','Freq_lo','Int_lo'};

for ss = 1:length(Participant)
%     try
%         Participant{ss}.name = Participant{ss}.namepostmerge;
%     end
    megpath{ss} = [pathstem Participant{ss}.groupfolder '/' Participant{ss}.namepostmerge '/' prefix Participant{ss}.namepostmerge '.mat'];
    if ~exist(megpath{ss},'file')
        megpath{ss} = [pathstem Participant{ss}.groupfolder '/' Participant{ss}.namepostmerge '/' prefix Participant{ss}.name{1} '.mat'];
    end
    diagnosis{ss} = Participant{ss}.diag;
   
end

[groups,~, group_inds] = unique(diagnosis,'stable');

if ~isempty(varargin)
    disp(['All groups in dataset'])
    groups
    [groups, group_ids, ~]=intersect(groups,varargin{1},'stable');
    disp(['Groups trimmed to'])
    groups
end

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
p_thresh = 0.05;

for j = 1:length(conditions)
    addpath('./stdshade')
    
    
    DEV_plot = figure(20*j);
    
        
        hold on
        for grp = 1:length(groups)
            linehandle(grp) = stdshade_TEC(squeeze(all_DEV(group_inds==group_ids(grp),:,j)),0.2,cmap(grp,:),D{1}.time,1,1);
        end
        % Plot statistically significant difference from control (FDR corrected for multiple comparisons across time)
        these_ylims = ylim;
        this_ylim_range = these_ylims(2)-these_ylims(1);
        for grp = 2:length(groups)
            for t = 1:size(D{1}.time,2)
                [h_bytime(t),p_bytime(t)]=ttest2(squeeze(all_DEV(group_inds==1,t,j)),squeeze(all_DEV(group_inds==group_ids(grp),t,j)));
            end
            all_fdr_p = mafdr(p_bytime(D{1}.time>0));
            all_positive_times = D{1}.time(D{1}.time>0);
            if any(all_positive_times(all_fdr_p<p_thresh))
                plot(all_positive_times(all_fdr_p<p_thresh),these_ylims(2)+(grp*this_ylim_range/40),'.','Color',cmap(grp,:))
            end
        end
        ylim([these_ylims(1),these_ylims(2)+(this_ylim_range/8)])
        
        xlabel('Time (s)')
        xlim([-0.1 0.500])
        plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
        ylabel([conditions{j} ' Response (AU)'])
        set(gca,'FontWeight','bold')
        set(gca,'LineWidth',2)
        
            legend(linehandle,groups)
        saveas(DEV_plot,['./outputfigures/' conditions{j} ' Response (AU).png']);
    
    
    Abs_DEV_plot = figure(400*j);
    
        
        hold on
        for grp = 1:length(groups)
            linehandle(grp) = stdshade_TEC(squeeze(abs(all_DEV(group_inds==group_ids(grp),:,j))),0.2,cmap(grp,:),D{1}.time,1,1);
        end
        % Plot statistically significant difference from control (FDR corrected for multiple comparisons across time)
        these_ylims = ylim;
        this_ylim_range = these_ylims(2)-these_ylims(1);
        for grp = 2:length(groups)
            for t = 1:size(D{1}.time,2)
                [h_bytime(t),p_bytime(t)]=ttest2(squeeze(abs(all_DEV(group_inds==1,t,j))),squeeze(abs(all_DEV(group_inds==group_ids(grp),t,j))));
            end
            all_fdr_p = mafdr(p_bytime(D{1}.time>0));
            all_positive_times = D{1}.time(D{1}.time>0);
            if any(all_positive_times(all_fdr_p<p_thresh))
                plot(all_positive_times(all_fdr_p<p_thresh),these_ylims(2)+(grp*this_ylim_range/40),'.','Color',cmap(grp,:))
            end
        end
        ylim([these_ylims(1),these_ylims(2)+(this_ylim_range/8)])
        
        xlabel('Time (s)')
        xlim([-0.1 0.500])
        plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
        ylabel([conditions{j} ' Abs Response (AU)'])
        set(gca,'FontWeight','bold')
        set(gca,'LineWidth',2)
        
            legend(linehandle,groups)
        saveas(Abs_DEV_plot,['./outputfigures/' conditions{j} ' Abs Response (AU).png']);
    
    
    if j > 1 %Don't do MMN for STDs only
        MMN_plot = figure(8000*j);
        
            
            hold on
            for grp = 1:length(groups)
                linehandle(grp) = stdshade_TEC(squeeze(all_MMN(group_inds==group_ids(grp),:,j)),0.2,cmap(grp,:),D{1}.time,1,1);
            end
            % Plot statistically significant difference from control (FDR corrected for multiple comparisons across time)
            these_ylims = ylim;
            this_ylim_range = these_ylims(2)-these_ylims(1);
            for grp = 2:length(groups)
                for t = 1:size(D{1}.time,2)
                    [h_bytime(t),p_bytime(t)]=ttest2(squeeze(all_MMN(group_inds==1,t,j)),squeeze(all_MMN(group_inds==group_ids(grp),t,j)));
                end
                all_fdr_p = mafdr(p_bytime(D{1}.time>0));
                all_positive_times = D{1}.time(D{1}.time>0);
                if any(all_positive_times(all_fdr_p<p_thresh))
                    plot(all_positive_times(all_fdr_p<p_thresh),these_ylims(2)+(grp*this_ylim_range/40),'.','Color',cmap(grp,:))
                end
            end
            ylim([these_ylims(1),these_ylims(2)+(this_ylim_range/8)])
            
            xlabel('Time (s)')
            xlim([-0.1 0.500])
            plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
            ylabel([conditions{j} ' Mismatch Response (AU)'])
            set(gca,'FontWeight','bold')
            set(gca,'LineWidth',2)
            
                legend(linehandle,groups)
            saveas(MMN_plot,['./outputfigures/' conditions{j} ' Mismatch Response (AU).png']);
        
        
        ABS_MMN_plot = figure(160000*j);
        
            
            hold on
            for grp = 1:length(groups)
                linehandle(grp) = stdshade_TEC(squeeze(all_abs_MMN(group_inds==group_ids(grp),:,j)),0.2,cmap(grp,:),D{1}.time,1,1);
            end
            % Plot statistically significant difference from control (FDR corrected for multiple comparisons across time)
            these_ylims = ylim;
            this_ylim_range = these_ylims(2)-these_ylims(1);
            for grp = 2:length(groups)
                for t = 1:size(D{1}.time,2)
                    [h_bytime(t),p_bytime(t)]=ttest2(squeeze(all_abs_MMN(group_inds==1,t,j)),squeeze(all_abs_MMN(group_inds==group_ids(grp),t,j)));
                end
                all_fdr_p = mafdr(p_bytime(D{1}.time>0));
                all_positive_times = D{1}.time(D{1}.time>0);
                if any(all_positive_times(all_fdr_p<p_thresh))
                    plot(all_positive_times(all_fdr_p<p_thresh),these_ylims(2)+(grp*this_ylim_range/40),'.','Color',cmap(grp,:))
                end
            end
            ylim([these_ylims(1),these_ylims(2)+(this_ylim_range/8)])
            
            xlabel('Time (s)')
            xlim([-0.1 0.500])
            plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
            ylabel([conditions{j} ' Abs Mismatch Response (AU)'])
            set(gca,'FontWeight','bold')
            set(gca,'LineWidth',2)
            
                legend(linehandle,groups)
            saveas(ABS_MMN_plot,['./outputfigures/' conditions{j} ' Abs Mismatch Response (AU).png']);
        
    end
    %pause
    close all %To prevent Java memory error
    
%     DEV_plot = figure(20*j);
% 
%         hold on
%         for grp = 1:length(groups)
%             linehandle(grp) = stdshade_TEC(squeeze(all_DEV(group_inds==group_ids(grp),:,j)),0.2,cmap(grp,:),D{1}.time,1,1);
%         end
%         xlabel('Time (s)')
%         xlim([-0.1 0.500])
%         plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
%         ylabel([conditions{j} ' Response (AU)'])
%         set(gca,'FontWeight','bold')
%         set(gca,'LineWidth',2)
%             legend(linehandle,groups)
%             
%             
%     
%     MMN_plot = figure(j);
% 
%         hold on
%         for grp = 1:length(groups)
%             linehandle(grp) = stdshade_TEC(squeeze(all_MMN(group_inds==group_ids(grp),:,j)),0.2,cmap(grp,:),D{1}.time,1,1);
%         end
%         xlabel('Time (s)')
%         xlim([-0.1 0.500])
%         plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
%         ylabel([conditions{j} ' Mismatch Response (AU)'])
%         set(gca,'FontWeight','bold')
%         set(gca,'LineWidth',2)
%             legend(linehandle,groups)
%     
%     ABS_MMN_plot = figure(400*j);
%         hold on
%         for grp = 1:length(groups)
%             linehandle(grp) = stdshade_TEC(squeeze(all_abs_MMN(group_inds==group_ids(grp),:,j)),0.2,cmap(grp,:),D{1}.time,1,1);
%         end
%         xlabel('Time (s)')
%         xlim([-0.1 0.500])
%         plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
%         ylabel([conditions{j} ' Abs Mismatch Response (AU)'])
%         set(gca,'FontWeight','bold')
%         set(gca,'LineWidth',2)
%             legend(linehandle,groups)
%             
%             DEV_indiv_plot = figure(8000*j);
% 
%         hold on
%         for grp = 1:length(groups)
%             linehandles{grp} = plot(D{1}.time,squeeze(all_DEV(group_inds==group_ids(grp),:,j)),'Color',cmap(grp,:));
%         end
%         xlabel('Time (s)')
%         xlim([-0.1 0.500])
%         ylabel([conditions{j} ' Individual Response (AU)'])
%         set(gca,'FontWeight','bold')
%         set(gca,'LineWidth',2)
%         for i = 1:length(linehandles)
%             linehandle(i) = linehandles{i}(1);
%         end
%             legend(linehandle,groups)
%     pause
%     close all %To prevent Java memory error
end
%pause
