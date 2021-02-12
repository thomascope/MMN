function quantify_MMN_LFP(Participant,pathstem,p,prefix,baselined)
% A function for plotting my extracted LFPs then quantifying various
% parameters after Hughes et al. 2013, namely:
% 1) Peak negative amplitude within 50-150 msec for the standard tone waveform
% 2) The mean amplitude of the MMN - originally 100-200ms but my data seems
% slightly delayed compared to Laura's, so going for 50ms either side of
% the first peak
% 3) First peak latency of the MMN 100-300ms (changed from 200ms to account
% for one of two quite delayed people and an overall right shift in data)
% 4) Ratio of 2:1 - to account for arbitrary units and relatively scale the MMN
addpath(['/group/language/data/thomascope/MMN/ICA_denoise/stdshade']);

Sname = {'left A1';
    'left STG';
    'left IFG';
    'left IPC';
    'right A1';
    'right STG';
    'right IFG';
    'right IPC'};

conditions = {'STD','DVT','Loc','Int','Dur','Gap','Freq'};

mkdir('./outputfigures/source_flipped/stats/')

if ~exist('baselined', 'var')
    baselined = 0;
end

for ss = 1:length(Participant)
    try
       % Participant{ss}.name = Participant{ss}.namepostmerge;
    end
    
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
    all_times{ss}= D{ss}.time;
    for i = 1:length(Sname)
        %First check if data needs to be flipped - make all STD M100s in A1,
        %STG and IPC negative, and all second deflections in IFG positive
        %(M100 not always reliably seen in IFG)
        all_STD(i,ss,:)=D{ss}(i,:,1);
        if i == 3 || i == 7
            if abs(max(all_STD(i,ss,all_times{ss}>=0.09&all_times{ss}<=0.165))) < abs(min(all_STD(i,ss,all_times{ss}>=0.09&all_times{ss}<=0.165)))
                D{ss}(i,:,:) = -D{ss}(i,:,:);
                all_STD(i,ss,:)=D{ss}(i,:,1);
            end
        else
            if abs(max(all_STD(i,ss,all_times{ss}>=0.03&all_times{ss}<=0.08))) > abs(min(all_STD(i,ss,all_times{ss}>=0.03&all_times{ss}<=0.08)))
                D{ss}(i,:,:) = -D{ss}(i,:,:);
                all_STD(i,ss,:)=D{ss}(i,:,1);
            end
        end

        for j = 1:length(conditions)
            all_DEV(i,ss,:,j)=D{ss}(i,:,j);
            all_MMN(i,ss,:,j)=D{ss}(i,:,1)-D{ss}(i,:,j);
            all_abs_MMN(i,ss,:,j)=abs(all_MMN(i,ss,:,j));
        end
    end
% %     %Correct for changing sound system latency over the years - 13ms until
% %     %2011, 32ms until 2015, then 26ms - no longer needed, integrated
% %     % pipeline accounts for this
% %     temp_split = strsplit(Participant{ss}.name,'meg');
% %     this_year(ss) = str2num(temp_split{2}(1:2));
% %     if this_year(ss)<=11
% %         all_times{ss}= D{ss}.time-0.013;
% %     elseif this_year(ss)>=15
% %         all_times{ss} = D{ss}.time-0.026;
% %     else
% %         all_times{ss} = D{ss}.time-0.032;
% %     end
end


plot_points = 0;
for ss = 1:length(Participant)
    for i = 1:length(Sname)
        STD_M100_amplitude(i,ss)=max(abs(all_STD(i,ss,all_times{ss}>=0.05&all_times{ss}<=0.15))); %STD amplitude
        if plot_points == 1
            values_plot = figure(1);
            clf(gcf)
            set(gcf,'Position',[100 100 1600 800]);
            set(gcf, 'PaperPositionMode', 'auto');
            hold off
            subplot(2,ceil(length(conditions)/2),1)
            plot(all_times{ss},squeeze(all_STD(i,ss,:)),'k-')
            hold on
            plot(all_times{ss},repmat(STD_M100_amplitude(i,ss),1,length(all_times{ss})),'r--');
            plot(all_times{ss},-repmat(STD_M100_amplitude(i,ss),1,length(all_times{ss})),'r--');
            hold off
        end
        for j = 2:length(conditions)
            [~,firstpeakloc]=findpeaks(squeeze(all_abs_MMN(i,ss,all_times{ss}>=0.1&all_times{ss}<=0.25,j)),'MinPeakHeight',max(squeeze(all_abs_MMN(i,ss,all_times{ss}>=0.1&all_times{ss}<=0.25,j)))*3/4);
            if isempty(firstpeakloc)
                [~,firstpeakloc]=findpeaks(squeeze(all_abs_MMN(i,ss,all_times{ss}>=0.1&all_times{ss}<=0.25,j)),'MinPeakHeight',max(squeeze(all_abs_MMN(i,ss,all_times{ss}>=0.1&all_times{ss}<=0.25,j)))*0.5);
                if isempty(firstpeakloc)
                    [~,firstpeakloc] = max(squeeze(all_abs_MMN(i,ss,all_times{ss}>=0.1&all_times{ss}<=0.25,j)));
                end
            end
            firstpeaksample = min(find(all_times{ss}>=0.1))+firstpeakloc(1)-1;
            MMN_latency(i,ss,j)=all_times{ss}(firstpeaksample);
            MMN_amplitude(i,ss,j)=mean(all_abs_MMN(i,ss,all_times{ss}>=0.1&all_times{ss}<=0.2,j),3);
            %relative_MMN_amplitude(i,ss,j)=MMN_amplitude(i,ss,j)/STD_M100_amplitude(i,ss);
            relative_MMN_amplitude(i,ss,j)=mean(all_abs_MMN(i,ss,all_times{ss}>=(MMN_latency(i,ss,j)-0.05)&all_times{ss}<=(MMN_latency(i,ss,j)+0.05),j),3);
            if plot_points == 1
               subplot(2,ceil(length(conditions)/2),j) 
               plot(all_times{ss},squeeze(all_abs_MMN(i,ss,:,j)),'k-')
               hold on
               plot(all_times{ss},squeeze(all_MMN(i,ss,:,j)),'k--')
               plot(MMN_latency(i,ss,j),MMN_amplitude(i,ss,j),'x')
               plot(MMN_latency(i,ss,j),-MMN_amplitude(i,ss,j),'x')
               drawnow
               hold off
            end
        end
        %pause
    end
end

cmap = colormap(parula(length(groups)));
p_thresh = 0.05;

M100_plot = figure(20000*j);
set(gcf,'Position',[100 100 1600 800]);
set(gcf, 'PaperPositionMode', 'auto');
M100_table = [table(diagnosis','VariableNames',{'Diagnosis'}),array2table(STD_M100_amplitude')];
M100_covariates = table(Sname,'VariableNames',{'Region'});
M100_rm = fitrm(M100_table,'Var1-Var8~Diagnosis','WithinDesign',M100_covariates);
M100_ranovatbl = ranova(M100_rm,'WithinModel','Region')
supertitle = [];
if M100_ranovatbl.pValueGG(find(strcmp(M100_ranovatbl.Row,'Diagnosis'))) < p_thresh
    supertitle = [supertitle 'Main effect of Diagnosis p=' num2str(M100_ranovatbl.pValueGG(find(strcmp(M100_ranovatbl.Row,'Diagnosis'))),3) '. '];
    Diagnosis_table = multcompare(M100_rm,'Diagnosis');
    significant_contrasts = Diagnosis_table.pValue<0.05;
    Diagnosis_table(significant_contrasts,:)
end
if M100_ranovatbl.pValueGG(find(strcmp(M100_ranovatbl.Row,'Diagnosis:Region'))) < p_thresh
    supertitle = [supertitle 'Diagnosis by region interaction p=' num2str(M100_ranovatbl.pValueGG(find(strcmp(M100_ranovatbl.Row,'Diagnosis:Region'))),3) '. '];
end
if isempty(supertitle)
    supertitle = [supertitle 'No main effect of Diagnosis or Diagnosis by region interaction in STD amplitude'];
end
for i = 1:8
    subplot(2,4,i)
    sgtitle(supertitle)
    hold on
    for grp = 1:length(groups)
        scatter(repmat(grp+0.1,1,sum(group_inds==grp)),STD_M100_amplitude(i,group_inds==grp))

        [~,p] = ttest2(STD_M100_amplitude(i,group_inds==1),STD_M100_amplitude(i,group_inds==grp));
        %if p<(p_thresh/(length(groups)-1)) %Bonferroni correct 2 sample t-tests vs controls
        if M100_ranovatbl.pValueGG(find(strcmp(M100_ranovatbl.Row,'Diagnosis:Region'))) < p_thresh && p<p_thresh %No longer Bonferroni as only do if diagnosis by region interaction is significant
            errorbar(grp-0.1,mean(STD_M100_amplitude(i,group_inds==grp),2),std(STD_M100_amplitude(i,group_inds==grp))/sqrt(sum(group_inds==grp)),'-s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1,'color','black')
        else
            errorbar(grp-0.1,mean(STD_M100_amplitude(i,group_inds==grp),2),std(STD_M100_amplitude(i,group_inds==grp))/sqrt(sum(group_inds==grp)),'-s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1,'color','black')
        end
    end
    xlim([0 length(groups)+1])
    xticks([1:length(groups)])
    xticklabels(groups)
    xtickangle(30)
    title(Sname{i},'FontSize',34)
    ylabel('STD Amplitude (AU)')
end
saveas(M100_plot,['./outputfigures/source_flipped/stats/M100 Amplitude (AU)_nolegend.png']);

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
            all_fdr_p = mafdr(p_bytime(D{1}.time>0),'BHFDR',true);
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
            %legend(linehandle,groups)
        end
    end
    saveas(DEV_plot,['./outputfigures/source_flipped/' conditions{j} ' Response (AU)_nolegend.png']);
    
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
            all_fdr_p = mafdr(p_bytime(D{1}.time>0),'BHFDR',true);
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
            %legend(linehandle,groups)
        end
    end
    saveas(Abs_DEV_plot,['./outputfigures/source_flipped/' conditions{j} ' Abs Response (AU)_nolegend.png']);
    
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
                all_fdr_p = mafdr(p_bytime(D{1}.time>0),'BHFDR',true);
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
                %legend(linehandle,groups)
            end
        end
        saveas(MMN_plot,['./outputfigures/source_flipped/' conditions{j} ' Mismatch Response (AU)_nolegend.png']);
        
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
                all_fdr_p = mafdr(p_bytime(D{1}.time>0),'BHFDR',true);
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
                %legend(linehandle,groups)
            end
        end
        saveas(ABS_MMN_plot,['./outputfigures/source_flipped/' conditions{j} ' Abs Mismatch Response (AU)_nolegend.png']);
    end
    %pause
%     close all %To prevent Java memory error
end
%Now plot each region by condition
for i = 1:8
    addpath('./stdshade')
    
    MMN_latency_plot = figure(20000*i);
    set(gcf,'Position',[100 100 1600 800]);
    set(gcf, 'PaperPositionMode', 'auto');
    
    latency_table = [table(diagnosis','VariableNames',{'Diagnosis'}),array2table(squeeze(MMN_latency(i,:,3:7)))];
    latency_covariates = table(conditions(3:7)','VariableNames',{'Condition'});
    latency_rm = fitrm(latency_table,'Var1-Var5~Diagnosis','WithinDesign',latency_covariates);
    latency_ranovatbl = ranova(latency_rm,'WithinModel','Condition')
    supertitle = [];
    if latency_ranovatbl.pValueGG(find(strcmp(latency_ranovatbl.Row,'Diagnosis'))) < p_thresh
        supertitle = [supertitle 'Latency main effect of Diagnosis p=' num2str(latency_ranovatbl.pValueGG(find(strcmp(latency_ranovatbl.Row,'Diagnosis'))),3) ' in ' Sname{i} '. '];
        Diagnosis_table = multcompare(latency_rm,'Diagnosis');
        significant_contrasts = Diagnosis_table.pValue<0.05;
        Sname{i}
        Diagnosis_table(significant_contrasts,:)
    end
    if latency_ranovatbl.pValueGG(find(strcmp(latency_ranovatbl.Row,'Diagnosis:Condition'))) < p_thresh
        supertitle = [supertitle 'Latency diagnosis by condition interaction p=' num2str(latency_ranovatbl.pValueGG(find(strcmp(latency_ranovatbl.Row,'Diagnosis:Condition'))),3) ' in ' Sname{i} '. '];
    end
    if anova1(squeeze(MMN_latency(i,:,2)),diagnosis,'off') < p_thresh
        supertitle = [supertitle 'One way ANOVA group difference for DVT MMN latency in ' Sname{i} ' p=' num2str(anova1(squeeze(relative_MMN_amplitude(i,:,2)),diagnosis,'off'),3) '. '];
    end
    if isempty(supertitle)
        supertitle = [supertitle 'No main effect of Diagnosis or Diagnosis by region interaction in MMN latency in ' Sname{i}];
    end
    
    for j = 2:7
        subplot(2,3,j-1)
        sgtitle(supertitle)
        hold on
        for grp = 1:length(groups)
            scatter(repmat(grp+0.1,1,sum(group_inds==grp)),squeeze(MMN_latency(i,group_inds==grp,j)));
            [~,p] = ttest2(squeeze(MMN_latency(i,group_inds==1,j)),squeeze(MMN_latency(i,group_inds==grp,j)));
            %if p<(p_thresh/(length(groups)-1)) %Bonferroni correct 2 sample t-tests vs controls
            if j == 2
                if anova1(squeeze(MMN_latency(i,:,j)),diagnosis,'off') < p_thresh && p<p_thresh %No longer Bonferroni as only do if diagnosis by region interaction is significant
                    errorbar(grp-0.1,mean(squeeze(MMN_latency(i,group_inds==grp,j)),2),std(squeeze(MMN_latency(i,group_inds==grp,j)))/sqrt(sum(group_inds==grp)),'-s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1,'color','black')
                else
                    errorbar(grp-0.1,mean(squeeze(MMN_latency(i,group_inds==grp,j)),2),std(squeeze(MMN_latency(i,group_inds==grp,j)))/sqrt(sum(group_inds==grp)),'-s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1,'color','black')
                end
            else
                if latency_ranovatbl.pValueGG(find(strcmp(latency_ranovatbl.Row,'Diagnosis:Condition'))) < p_thresh && p<p_thresh %No longer Bonferroni as only do if diagnosis by region interaction is significant
                    errorbar(grp-0.1,mean(squeeze(MMN_latency(i,group_inds==grp,j)),2),std(squeeze(MMN_latency(i,group_inds==grp,j)))/sqrt(sum(group_inds==grp)),'-s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1,'color','black')
                else
                    errorbar(grp-0.1,mean(squeeze(MMN_latency(i,group_inds==grp,j)),2),std(squeeze(MMN_latency(i,group_inds==grp,j)))/sqrt(sum(group_inds==grp)),'-s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1,'color','black')
                end
            end
        end
        xlim([0 length(groups)+1])
        xticks([1:length(groups)])
        xticklabels(groups)
        xtickangle(30)
        title(conditions{j},'FontSize',34)
        ylabel('MMN Latency (s)')
    end
    saveas(MMN_latency_plot,['./outputfigures/source_flipped/stats/' Sname{i} ' MMN Latency (s)_nolegend.png']);
    
        
    relative_MMN_amplitude_plot = figure(2000000*i);
    set(gcf,'Position',[100 100 1600 800]);
    set(gcf, 'PaperPositionMode', 'auto');
    
    amplitude_table = [table(diagnosis','VariableNames',{'Diagnosis'}),array2table(squeeze(relative_MMN_amplitude(i,:,3:7)))];
    amplitude_covariates = table(conditions(3:7)','VariableNames',{'Condition'});
    amplitude_rm = fitrm(amplitude_table,'Var1-Var5~Diagnosis','WithinDesign',amplitude_covariates);
    relative_amplitude_ranovatbl = ranova(amplitude_rm,'WithinModel','Condition')
    supertitle = [];
    if relative_amplitude_ranovatbl.pValueGG(find(strcmp(relative_amplitude_ranovatbl.Row,'Diagnosis'))) < p_thresh
        supertitle = [supertitle 'Amplitude main effect of Diagnosis p=' num2str(relative_amplitude_ranovatbl.pValueGG(find(strcmp(relative_amplitude_ranovatbl.Row,'Diagnosis'))),3) ' in ' Sname{i} '. '];
        Diagnosis_table = multcompare(amplitude_rm,'Diagnosis');
        significant_contrasts = Diagnosis_table.pValue<0.05;
        Sname{i}
        Diagnosis_table(significant_contrasts,:)
    end
    if relative_amplitude_ranovatbl.pValueGG(find(strcmp(relative_amplitude_ranovatbl.Row,'Diagnosis:Condition'))) < p_thresh
        supertitle = [supertitle 'Amplitude diagnosis by condition interaction p=' num2str(relative_amplitude_ranovatbl.pValueGG(find(strcmp(relative_amplitude_ranovatbl.Row,'Diagnosis:Condition'))),3) ' in ' Sname{i} '. '];
    end
    if anova1(squeeze(relative_MMN_amplitude(i,:,2)),diagnosis,'off') < p_thresh
        supertitle = [supertitle 'One way ANOVA group difference for DVT MMN amplitude in ' Sname{i} ' p=' num2str(anova1(squeeze(relative_MMN_amplitude(i,:,2)),diagnosis,'off'),3) '. '];
    end
    if isempty(supertitle)
        supertitle = [supertitle 'No main effect of Diagnosis or Diagnosis by region interaction in MMN amplitude in ' Sname{i}];
    end
    
    for j = 2:7
        subplot(2,3,j-1)
        hold on
        sgtitle(supertitle)
        for grp = 1:length(groups)
            scatter(repmat(grp+0.1,1,sum(group_inds==grp)),squeeze(relative_MMN_amplitude(i,group_inds==grp,j)));
            [~,p] = ttest2(squeeze(relative_MMN_amplitude(i,group_inds==1,j)),squeeze(relative_MMN_amplitude(i,group_inds==grp,j)));
            %if p<(p_thresh/(length(groups)-1)) %Bonferroni correct 2 sample t-tests vs controls
            if j == 2
                if anova1(squeeze(relative_MMN_amplitude(i,:,j)),diagnosis,'off') < p_thresh && p<p_thresh %No longer Bonferroni as only do if diagnosis by region interaction is significant
                    errorbar(grp-0.1,mean(squeeze(relative_MMN_amplitude(i,group_inds==grp,j)),2),std(squeeze(relative_MMN_amplitude(i,group_inds==grp,j)))/sqrt(sum(group_inds==grp)),'-s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1,'color','black')
                else
                    errorbar(grp-0.1,mean(squeeze(relative_MMN_amplitude(i,group_inds==grp,j)),2),std(squeeze(relative_MMN_amplitude(i,group_inds==grp,j)))/sqrt(sum(group_inds==grp)),'-s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1,'color','black')
                end
            else
                if relative_amplitude_ranovatbl.pValueGG(find(strcmp(relative_amplitude_ranovatbl.Row,'Diagnosis:Condition'))) < p_thresh && p<p_thresh %No longer Bonferroni as only do if diagnosis by region interaction is significant
                    errorbar(grp-0.1,mean(squeeze(relative_MMN_amplitude(i,group_inds==grp,j)),2),std(squeeze(relative_MMN_amplitude(i,group_inds==grp,j)))/sqrt(sum(group_inds==grp)),'-s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1,'color','black')
                else
                    errorbar(grp-0.1,mean(squeeze(relative_MMN_amplitude(i,group_inds==grp,j)),2),std(squeeze(relative_MMN_amplitude(i,group_inds==grp,j)))/sqrt(sum(group_inds==grp)),'-s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1,'color','black')
                end
            end
        end
        xlim([0 length(groups)+1])
        xticks([1:length(groups)])
        xticklabels(groups)
        xtickangle(30)
        title(conditions{j},'FontSize',34)
        ylabel('Relative MMN amplitude')
    end
    saveas(relative_MMN_amplitude_plot,['./outputfigures/source_flipped/stats/' Sname{i} ' relative MMN amplitude_nolegend.png']);
    
    MMN_amplitude_plot = figure(2000001*i);
    set(gcf,'Position',[100 100 1600 800]);
    set(gcf, 'PaperPositionMode', 'auto');
    
    amplitude_table = [table(diagnosis','VariableNames',{'Diagnosis'}),array2table(squeeze(MMN_amplitude(i,:,3:7)))];
    amplitude_covariates = table(conditions(3:7)','VariableNames',{'Condition'});
    amplitude_rm = fitrm(amplitude_table,'Var1-Var5~Diagnosis','WithinDesign',amplitude_covariates);
    amplitude_ranovatbl = ranova(amplitude_rm,'WithinModel','Condition')
    supertitle = [];
    if amplitude_ranovatbl.pValueGG(find(strcmp(amplitude_ranovatbl.Row,'Diagnosis'))) < p_thresh
        supertitle = [supertitle 'Amplitude main effect of Diagnosis p=' num2str(amplitude_ranovatbl.pValueGG(find(strcmp(amplitude_ranovatbl.Row,'Diagnosis'))),3) ' in ' Sname{i} '. '];
        Diagnosis_table = multcompare(amplitude_rm,'Diagnosis');
        significant_contrasts = Diagnosis_table.pValue<0.05;
        Sname{i}
        Diagnosis_table(significant_contrasts,:)
    end
    if amplitude_ranovatbl.pValueGG(find(strcmp(amplitude_ranovatbl.Row,'Diagnosis:Condition'))) < p_thresh
        supertitle = [supertitle 'Amplitude diagnosis by condition interaction p=' num2str(amplitude_ranovatbl.pValueGG(find(strcmp(amplitude_ranovatbl.Row,'Diagnosis:Condition'))),3) ' in ' Sname{i} '. '];
    end
    if anova1(squeeze(MMN_amplitude(i,:,2)),diagnosis,'off') < p_thresh
        supertitle = [supertitle 'One way ANOVA group difference for DVT MMN amplitude in ' Sname{i} ' p=' num2str(anova1(squeeze(MMN_amplitude(i,:,2)),diagnosis,'off'),3) '. '];
    end
    if isempty(supertitle)
        supertitle = [supertitle 'No main effect of Diagnosis or Diagnosis by region interaction in MMN amplitude in ' Sname{i}];
    end
    
    for j = 2:7
        subplot(2,3,j-1)
        hold on
        sgtitle(supertitle)
        for grp = 1:length(groups)
            scatter(repmat(grp+0.1,1,sum(group_inds==grp)),squeeze(MMN_amplitude(i,group_inds==grp,j)));
            [~,p] = ttest2(squeeze(MMN_amplitude(i,group_inds==1,j)),squeeze(MMN_amplitude(i,group_inds==grp,j)));
            %if p<(p_thresh/(length(groups)-1)) %Bonferroni correct 2 sample t-tests vs controls
            if j == 2
                if anova1(squeeze(MMN_amplitude(i,:,j)),diagnosis,'off') < p_thresh && p<p_thresh %No longer Bonferroni as only do if diagnosis by region interaction is significant
                    errorbar(grp-0.1,mean(squeeze(MMN_amplitude(i,group_inds==grp,j)),2),std(squeeze(MMN_amplitude(i,group_inds==grp,j)))/sqrt(sum(group_inds==grp)),'-s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1,'color','black')
                else
                    errorbar(grp-0.1,mean(squeeze(MMN_amplitude(i,group_inds==grp,j)),2),std(squeeze(MMN_amplitude(i,group_inds==grp,j)))/sqrt(sum(group_inds==grp)),'-s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1,'color','black')
                end
            else
                if amplitude_ranovatbl.pValueGG(find(strcmp(amplitude_ranovatbl.Row,'Diagnosis:Condition'))) < p_thresh && p<p_thresh %No longer Bonferroni as only do if diagnosis by region interaction is significant
                    errorbar(grp-0.1,mean(squeeze(MMN_amplitude(i,group_inds==grp,j)),2),std(squeeze(MMN_amplitude(i,group_inds==grp,j)))/sqrt(sum(group_inds==grp)),'-s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1,'color','black')
                else
                    errorbar(grp-0.1,mean(squeeze(MMN_amplitude(i,group_inds==grp,j)),2),std(squeeze(MMN_amplitude(i,group_inds==grp,j)))/sqrt(sum(group_inds==grp)),'-s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1,'color','black')
                end
            end
        end
        xlim([0 length(groups)+1])
        xticks([1:length(groups)])
        xticklabels(groups)
        xtickangle(30)
        title(conditions{j},'FontSize',34)
        ylabel('MMN amplitude (AU)')
    end
    saveas(MMN_amplitude_plot,['./outputfigures/source_flipped/stats/' Sname{i} ' MMN amplitude (AU)_nolegend.png']);
    
    all_pvals_maineffect(i,:) = [latency_ranovatbl.pValueGG(find(strcmp(latency_ranovatbl.Row,'Diagnosis'))), relative_amplitude_ranovatbl.pValueGG(find(strcmp(relative_amplitude_ranovatbl.Row,'Diagnosis'))), amplitude_ranovatbl.pValueGG(find(strcmp(relative_amplitude_ranovatbl.Row,'Diagnosis')))];
    all_pvals_interaction(i,:) = [latency_ranovatbl.pValueGG(find(strcmp(latency_ranovatbl.Row,'Diagnosis:Condition'))), relative_amplitude_ranovatbl.pValueGG(find(strcmp(relative_amplitude_ranovatbl.Row,'Diagnosis:Condition'))), amplitude_ranovatbl.pValueGG(find(strcmp(relative_amplitude_ranovatbl.Row,'Diagnosis:Condition')))];
    all_pvals_key = {'latency','amplitude_around_peak','amplitude_100-200'};
    
%     DEV_plot = figure(20*i);
%     set(gcf,'Position',[100 100 1600 800]);
%     set(gcf, 'PaperPositionMode', 'auto');
%     for j = 1:7
%         subplot(2,4,j)
%         hold on
%         for grp = 1:length(groups)
%             linehandle(grp) = stdshade_TEC(squeeze(all_DEV(i,group_inds==grp,:,j)),0.2,cmap(grp,:),D{1}.time,1,1);
%         end
%         % Plot statistically significant difference from control (FDR corrected for multiple comparisons across time)
%         these_ylims = ylim;
%         this_ylim_range = these_ylims(2)-these_ylims(1);
%         for grp = 2:length(groups)
%             for t = 1:size(D{1}.time,2)
%                 [h_bytime(t),p_bytime(t)]=ttest2(squeeze(all_DEV(i,group_inds==1,t,j)),squeeze(all_DEV(i,group_inds==grp,t,j)));
%             end
%             all_fdr_p = mafdr(p_bytime(D{1}.time>0));
%             all_positive_times = D{1}.time(D{1}.time>0);
%             if any(all_positive_times(all_fdr_p<p_thresh))
%                 plot(all_positive_times(all_fdr_p<p_thresh),these_ylims(2)+(grp*this_ylim_range/40),'.','Color',cmap(grp,:))
%             end
%         end
%         ylim([these_ylims(1),these_ylims(2)+(this_ylim_range/8)])
%         title(conditions{j},'FontSize',34)
%         xlabel('Time (s)')
%         xlim([-0.1 0.500])
%         plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
%         ylabel([conditions{j} ' Response (AU)'])
%         set(gca,'FontWeight','bold')
%         set(gca,'LineWidth',2)
%     end
%     for j = 8
%         subplot(2,4,j)
%         hold on
%         for grp = 1:length(groups)
%             linehandle(grp) = stdshade_TEC(squeeze(all_DEV(i,group_inds==grp,:,1)),0,cmap(grp,:),D{1}.time,1,1);
%         end
%         axis off
%         legend(linehandle,groups)
%         drawnow
%         delete(linehandle)
%         drawnow
%     end
%     saveas(DEV_plot,['./outputfigures/source_flipped/' Sname{i} ' Response (AU)_nolegend.png']);
%     
%     Abs_DEV_plot = figure(400*i);
%     set(gcf,'Position',[100 100 1600 800]);
%     set(gcf, 'PaperPositionMode', 'auto');
%     for j = 1:7
%         subplot(2,4,j)
%         hold on
%         for grp = 1:length(groups)
%             linehandle(grp) = stdshade_TEC(squeeze(abs(all_DEV(i,group_inds==grp,:,j))),0.2,cmap(grp,:),D{1}.time,1,1);
%         end
%         % Plot statistically significant difference from control (FDR corrected for multiple comparisons across time)
%         these_ylims = ylim;
%         this_ylim_range = these_ylims(2)-these_ylims(1);
%         for grp = 2:length(groups)
%             for t = 1:size(D{1}.time,2)
%                 [h_bytime(t),p_bytime(t)]=ttest2(squeeze(abs(all_DEV(i,group_inds==1,t,j))),squeeze(abs(all_DEV(i,group_inds==grp,t,j))));
%             end
%             all_fdr_p = mafdr(p_bytime(D{1}.time>0));
%             all_positive_times = D{1}.time(D{1}.time>0);
%             if any(all_positive_times(all_fdr_p<p_thresh))
%                 plot(all_positive_times(all_fdr_p<p_thresh),these_ylims(2)+(grp*this_ylim_range/40),'.','Color',cmap(grp,:))
%             end
%         end
%         ylim([these_ylims(1),these_ylims(2)+(this_ylim_range/8)])
%         title(['Abs ' conditions{j}],'FontSize',28)
%         xlabel('Time (s)')
%         xlim([-0.1 0.500])
%         plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
%         ylabel([conditions{j} ' Abs Response (AU)'])
%         set(gca,'FontWeight','bold')
%         set(gca,'LineWidth',2)
%     end
%     for j = 8
%         subplot(2,4,j)
%         hold on
%         for grp = 1:length(groups)
%             linehandle(grp) = stdshade_TEC(squeeze(all_DEV(i,group_inds==grp,:,1)),0,cmap(grp,:),D{1}.time,1,1);
%         end
%         axis off
%         legend(linehandle,groups)
%         drawnow
%         delete(linehandle)
%         drawnow
%     end
%     saveas(Abs_DEV_plot,['./outputfigures/source_flipped/' Sname{i} ' Abs Response (AU)_nolegend.png']);
%     
%         MMN_plot = figure(8000*i);
%         set(gcf,'Position',[100 100 1600 800]);
%         set(gcf, 'PaperPositionMode', 'auto');
%         for j = 1:7
%             subplot(2,4,j)
%             hold on
%             for grp = 1:length(groups)
%                 linehandle(grp) = stdshade_TEC(squeeze(all_MMN(i,group_inds==grp,:,j)),0.2,cmap(grp,:),D{1}.time,1,1);
%             end
%             % Plot statistically significant difference from control (FDR corrected for multiple comparisons across time)
%             these_ylims = ylim;
%             this_ylim_range = these_ylims(2)-these_ylims(1);
%             if j >= 2
%             for grp = 2:length(groups)
%                 for t = 1:size(D{1}.time,2)
%                     [h_bytime(t),p_bytime(t)]=ttest2(squeeze(all_MMN(i,group_inds==1,t,j)),squeeze(all_MMN(i,group_inds==grp,t,j)));
%                 end
%                 all_fdr_p = mafdr(p_bytime(D{1}.time>0));
%                 all_positive_times = D{1}.time(D{1}.time>0);
%                 if any(all_positive_times(all_fdr_p<p_thresh))
%                     plot(all_positive_times(all_fdr_p<p_thresh),these_ylims(2)+(grp*this_ylim_range/40),'.','Color',cmap(grp,:))
%                 end
%             end
%             end
%             ylim([these_ylims(1),these_ylims(2)+(this_ylim_range/8)])
%             title(['MMN ' conditions{j}],'FontSize',28)
%             xlabel('Time (s)')
%             xlim([-0.1 0.500])
%             plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
%             ylabel([conditions{j} ' Mismatch Response (AU)'])
%             set(gca,'FontWeight','bold')
%             set(gca,'LineWidth',2)
%         end
%         for j = 8
%             subplot(2,4,j)
%             hold on
%             for grp = 1:length(groups)
%                 linehandle(grp) = stdshade_TEC(squeeze(all_MMN(i,group_inds==grp,:,1)),0,cmap(grp,:),D{1}.time,1,1);
%             end
%             axis off
%             legend(linehandle,groups)
%             drawnow
%             delete(linehandle)
%             drawnow
%         end
%         saveas(MMN_plot,['./outputfigures/source_flipped/' Sname{i} ' Mismatch Response (AU)_nolegend.png']);
%         
%         ABS_MMN_plot = figure(160000*i);
%         set(gcf,'Position',[100 100 1600 800]);
%         set(gcf, 'PaperPositionMode', 'auto');
%         for j = 1:7
%             subplot(2,4,j)
%             hold on
%             for grp = 1:length(groups)
%                 linehandle(grp) = stdshade_TEC(squeeze(all_abs_MMN(i,group_inds==grp,:,j)),0.2,cmap(grp,:),D{1}.time,1,1);
%             end
%             % Plot statistically significant difference from control (FDR corrected for multiple comparisons across time)
%             these_ylims = ylim;
%             this_ylim_range = these_ylims(2)-these_ylims(1);
%             if j >= 2
%             for grp = 2:length(groups)
%                 for t = 1:size(D{1}.time,2)
%                     [h_bytime(t),p_bytime(t)]=ttest2(squeeze(all_abs_MMN(i,group_inds==1,t,j)),squeeze(all_abs_MMN(i,group_inds==grp,t,j)));
%                 end
%                 all_fdr_p = mafdr(p_bytime(D{1}.time>0));
%                 all_positive_times = D{1}.time(D{1}.time>0);
%                 if any(all_positive_times(all_fdr_p<p_thresh))
%                     plot(all_positive_times(all_fdr_p<p_thresh),these_ylims(2)+(grp*this_ylim_range/40),'.','Color',cmap(grp,:))
%                 end
%             end
%             end
%             ylim([these_ylims(1),these_ylims(2)+(this_ylim_range/8)])
%             title(['Abs MMN ' conditions{j}],'FontSize',22)
%             xlabel('Time (s)')
%             xlim([-0.1 0.500])
%             plot(D{1}.time,zeros(1,length(D{1}.time)),'k--','LineWidth',2)
%             ylabel([conditions{j} ' Abs Mismatch Response (AU)'])
%             set(gca,'FontWeight','bold')
%             set(gca,'LineWidth',2)
%         end
%         for j = 8
%             subplot(2,4,j)
%             hold on
%             for grp = 1:length(groups)
%                 linehandle(grp) = stdshade_TEC(squeeze(all_abs_MMN(i,group_inds==grp,:,1)),0,cmap(grp,:),D{1}.time,1,1);
%             end
%             axis off
%             legend(linehandle,groups)
%             drawnow
%             delete(linehandle)
%             drawnow
%         end
%         saveas(ABS_MMN_plot,['./outputfigures/source_flipped/' Sname{i} ' Abs Mismatch Response (AU)_nolegend.png']);
% 
%     %pause
    close all %To prevent Java memory error
end
save(['./outputfigures/source_flipped/stats/all_lfp_pvals.mat'], 'all_pvals_maineffect', 'all_pvals_interaction', 'all_pvals_key')
    
pause
