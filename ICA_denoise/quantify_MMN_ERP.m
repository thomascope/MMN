function quantify_MMN_ERP(Participant,pathstem,p,prefix,varargin)
% A function for plotting my extracted ERPs then quantifying various
% parameters after Hughes et al. 2013, namely:
% 1) Peak negative amplitude within 50-150 msec for the standard tone waveform
% 2) The mean amplitude of the MMN - originally 100-200ms but my data seems
% slightly delayed compared to Laura's, so going for 50ms either side of
% the first peak
% 3) First peak latency of the MMN 100-300ms (changed from 200ms to account
% for one of two quite delayed people and an overall right shift in data)
% 4) Ratio of 2:1 - to account for arbitrary units and relatively scale the MMN
addpath(['/group/language/data/thomascope/MMN/ICA_denoise/stdshade']);

conditions = {'STD','DVT','Loc','Int','Dur','Gap','Freq'};

mkdir('./outputfigures/scalp/stats/')

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
        all_abs_MMN(ss,:,j)=abs(all_MMN(ss,:,j));
    end
end


plot_points = 1;
for ss = 1:length(Participant)
    STD_M100_amplitude(ss)=max(abs(all_STD(ss,D{ss}.time>=0.05&D{ss}.time<=0.15))); %STD amplitude
    if plot_points == 1
        values_plot = figure(1);
        clf(gcf)
        set(gcf,'Position',[100 100 1600 800]);
        set(gcf, 'PaperPositionMode', 'auto');
        hold off
        subplot(2,ceil(length(conditions)/2),1)
        plot(D{ss}.time,squeeze(all_STD(ss,:)),'k-')
        hold on
        plot(D{ss}.time,repmat(STD_M100_amplitude(ss),1,length(D{ss}.time)),'r--');
        plot(D{ss}.time,-repmat(STD_M100_amplitude(ss),1,length(D{ss}.time)),'r--');
        hold off
    end
    for j = 2:length(conditions)
        [~,firstpeakloc]=findpeaks(squeeze(all_abs_MMN(ss,D{ss}.time>=0.1&D{ss}.time<=0.3,j)),'MinPeakHeight',max(squeeze(all_abs_MMN(ss,D{ss}.time>=0.1&D{ss}.time<=0.3,j)))*3/4);
        if isempty(firstpeakloc)
            [~,firstpeakloc]=findpeaks(squeeze(all_abs_MMN(ss,D{ss}.time>=0.1&D{ss}.time<=0.3,j)),'MinPeakHeight',max(squeeze(all_abs_MMN(ss,D{ss}.time>=0.1&D{ss}.time<=0.3,j)))*0.5);
            if isempty(firstpeakloc)
                [~,firstpeakloc] = max(squeeze(all_abs_MMN(ss,D{ss}.time>=0.1&D{ss}.time<=0.3,j)));
            end
        end
        firstpeaksample = min(find(D{ss}.time>=0.1))+firstpeakloc(1)-1;
        MMN_latency(ss,j)=D{ss}.time(firstpeaksample);
        MMN_amplitude(ss,j)=mean(all_abs_MMN(ss,D{ss}.time>=(MMN_latency(ss,j)-0.05)&D{ss}.time<=(MMN_latency(ss,j)+0.05),j),2);
        relative_MMN_amplitude(ss,j)=MMN_amplitude(ss,j)/STD_M100_amplitude(ss);
        if plot_points == 1
            subplot(2,ceil(length(conditions)/2),j)
            plot(D{ss}.time,squeeze(all_abs_MMN(ss,:,j)),'k-')
            hold on
            plot(D{ss}.time,squeeze(all_MMN(ss,:,j)),'k--')
            plot(MMN_latency(ss,j),MMN_amplitude(ss,j),'x')
            plot(MMN_latency(ss,j),-MMN_amplitude(ss,j),'x')
            drawnow
            hold off
        end
    end
    %pause
end

cmap = colormap(parula(length(groups)));
p_thresh = 0.05;

M100_plot = figure(20000*j);
set(gcf,'Position',[100 100 1600 800]);
set(gcf, 'PaperPositionMode', 'auto');
supertitle = [];
if anova1(squeeze(STD_M100_amplitude),diagnosis,'off') < p_thresh
    supertitle = [supertitle 'Main effect of Diagnosis p=' num2str(anova1(squeeze(STD_M100_amplitude),diagnosis,'off'),3) '. '];
end
if isempty(supertitle)
    supertitle = [supertitle 'No main effect of Diagnosis in STD amplitude'];
end
sgtitle(supertitle)
hold on
for grp = 1:length(groups)
    scatter(repmat(grp+0.1,1,sum(group_inds==grp)),STD_M100_amplitude(group_inds==grp))
    
    [~,p] = ttest2(STD_M100_amplitude(group_inds==1),STD_M100_amplitude(group_inds==grp));
    %if p<(p_thresh/(length(groups)-1)) %Bonferroni correct 2 sample t-tests vs controls
    if STD_M100_amplitude < p_thresh & p<p_thresh %No longer Bonferroni as only do if diagnosis by region interaction is significant
        errorbar(grp-0.1,mean(STD_M100_amplitude(group_inds==grp)),std(STD_M100_amplitude(group_inds==grp))/sqrt(sum(group_inds==grp)),'-s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1,'color','black')
    else
        errorbar(grp-0.1,mean(STD_M100_amplitude(group_inds==grp)),std(STD_M100_amplitude(group_inds==grp))/sqrt(sum(group_inds==grp)),'-s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1,'color','black')
    end
end
xlim([0 length(groups)+1])
xticks([1:length(groups)])
xticklabels(groups)
xtickangle(30)
ylabel('Magnetometer Amplitude (fT/mm)')
saveas(M100_plot,['./outputfigures/scalp/stats/M100 Amplitude (AU).png']);

MMN_latency_plot = figure(20000);
set(gcf,'Position',[100 100 1600 800]);
set(gcf, 'PaperPositionMode', 'auto');

latency_table = [table(diagnosis','VariableNames',{'Diagnosis'}),array2table(squeeze(MMN_latency(:,3:7)))];
latency_covariates = table(conditions(3:7)','VariableNames',{'Condition'});
latency_rm = fitrm(latency_table,'Var1-Var5~Diagnosis','WithinDesign',latency_covariates);
latency_ranovatbl = ranova(latency_rm,'WithinModel','Condition')
supertitle = [];
if latency_ranovatbl.pValueGG(find(strcmp(latency_ranovatbl.Row,'Diagnosis'))) < p_thresh
    supertitle = [supertitle 'Latency main effect of Diagnosis p=' num2str(latency_ranovatbl.pValueGG(find(strcmp(latency_ranovatbl.Row,'Diagnosis'))),3) ' on scalp. '];
    Diagnosis_table = multcompare(latency_rm,'Diagnosis');
    significant_contrasts = Diagnosis_table.pValue<0.05;
    Sname{i}
    Diagnosis_table(significant_contrasts,:)
end
if latency_ranovatbl.pValueGG(find(strcmp(latency_ranovatbl.Row,'Diagnosis:Condition'))) < p_thresh
    supertitle = [supertitle 'Latency diagnosis by condition interaction p=' num2str(latency_ranovatbl.pValueGG(find(strcmp(latency_ranovatbl.Row,'Diagnosis:Condition'))),3) ' on scalp. '];
end
if anova1(squeeze(MMN_latency(:,2)),diagnosis,'off') < p_thresh
    supertitle = [supertitle 'One way ANOVA group difference for DVT MMN latency on scalp p=' num2str(anova1(squeeze(relative_MMN_amplitude(:,2)),diagnosis,'off'),3) '. '];
end
if isempty(supertitle)
    supertitle = [supertitle 'No main effect of Diagnosis or Diagnosis by region interaction in MMN latency in ' Sname{i}];
end

for j = 2:7
    subplot(2,3,j-1)
    sgtitle(supertitle)
    hold on
    for grp = 1:length(groups)
        scatter(repmat(grp+0.1,1,sum(group_inds==grp)),squeeze(MMN_latency(group_inds==grp,j)));
        [~,p] = ttest2(squeeze(MMN_latency(group_inds==1,j)),squeeze(MMN_latency(group_inds==grp,j)));
        %if p<(p_thresh/(length(groups)-1)) %Bonferroni correct 2 sample t-tests vs controls
        if j == 2
            if anova1(squeeze(MMN_latency(:,j)),diagnosis,'off') < p_thresh && p<p_thresh %No longer Bonferroni as only do if diagnosis by region interaction is significant
                errorbar(grp-0.1,mean(squeeze(MMN_latency(group_inds==grp,j)),1),std(squeeze(MMN_latency(group_inds==grp,j)))/sqrt(sum(group_inds==grp)),'-s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1,'color','black')
            else
                errorbar(grp-0.1,mean(squeeze(MMN_latency(group_inds==grp,j)),1),std(squeeze(MMN_latency(group_inds==grp,j)))/sqrt(sum(group_inds==grp)),'-s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1,'color','black')
            end
        else
            if latency_ranovatbl.pValueGG(find(strcmp(latency_ranovatbl.Row,'Diagnosis:Condition'))) < p_thresh && p<p_thresh %No longer Bonferroni as only do if diagnosis by region interaction is significant
                errorbar(grp-0.1,mean(squeeze(MMN_latency(group_inds==grp,j)),1),std(squeeze(MMN_latency(group_inds==grp,j)))/sqrt(sum(group_inds==grp)),'-s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1,'color','black')
            else
                errorbar(grp-0.1,mean(squeeze(MMN_latency(group_inds==grp,j)),1),std(squeeze(MMN_latency(group_inds==grp,j)))/sqrt(sum(group_inds==grp)),'-s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1,'color','black')
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
saveas(MMN_latency_plot,['./outputfigures/scalp/stats/MMN Latency (s).png']);


relative_MMN_amplitude_plot = figure(2000000);
set(gcf,'Position',[100 100 1600 800]);
set(gcf, 'PaperPositionMode', 'auto');

amplitude_table = [table(diagnosis','VariableNames',{'Diagnosis'}),array2table(squeeze(relative_MMN_amplitude(:,3:7)))];
amplitude_covariates = table(conditions(3:7)','VariableNames',{'Condition'});
amplitude_rm = fitrm(amplitude_table,'Var1-Var5~Diagnosis','WithinDesign',amplitude_covariates);
amplitude_ranovatbl = ranova(amplitude_rm,'WithinModel','Condition')
supertitle = [];
if amplitude_ranovatbl.pValueGG(find(strcmp(amplitude_ranovatbl.Row,'Diagnosis'))) < p_thresh
    supertitle = [supertitle 'Amplitude main effect of Diagnosis p=' num2str(amplitude_ranovatbl.pValueGG(find(strcmp(amplitude_ranovatbl.Row,'Diagnosis'))),3) ' on scalp. '];
    Diagnosis_table = multcompare(amplitude_rm,'Diagnosis');
    significant_contrasts = Diagnosis_table.pValue<0.05;
    Sname{i}
    Diagnosis_table(significant_contrasts,:)
end
if amplitude_ranovatbl.pValueGG(find(strcmp(amplitude_ranovatbl.Row,'Diagnosis:Condition'))) < p_thresh
    supertitle = [supertitle 'Amplitude diagnosis by condition interaction p=' num2str(amplitude_ranovatbl.pValueGG(find(strcmp(amplitude_ranovatbl.Row,'Diagnosis:Condition'))),3) ' on scalp. '];
end
if anova1(squeeze(relative_MMN_amplitude(:,2)),diagnosis,'off') < p_thresh
    supertitle = [supertitle 'One way ANOVA group difference for DVT MMN amplitude on scalp p=' num2str(anova1(squeeze(relative_MMN_amplitude(:,2)),diagnosis,'off'),3) '. '];
end
if isempty(supertitle)
    supertitle = [supertitle 'No main effect of Diagnosis or Diagnosis by region interaction in MMN amplitude on scalp'];
end

for j = 2:7
    subplot(2,3,j-1)
    hold on
    sgtitle(supertitle)
    for grp = 1:length(groups)
        scatter(repmat(grp+0.1,1,sum(group_inds==grp)),squeeze(relative_MMN_amplitude(group_inds==grp,j)));
        [~,p] = ttest2(squeeze(relative_MMN_amplitude(group_inds==1,j)),squeeze(relative_MMN_amplitude(group_inds==grp,j)));
        %if p<(p_thresh/(length(groups)-1)) %Bonferroni correct 2 sample t-tests vs controls
        if j == 2
            if anova1(squeeze(relative_MMN_amplitude(:,j)),diagnosis,'off') < p_thresh && p<p_thresh %No longer Bonferroni as only do if diagnosis by region interaction is significant
                errorbar(grp-0.1,mean(squeeze(relative_MMN_amplitude(group_inds==grp,j)),1),std(squeeze(relative_MMN_amplitude(group_inds==grp,j)))/sqrt(sum(group_inds==grp)),'-s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1,'color','black')
            else
                errorbar(grp-0.1,mean(squeeze(relative_MMN_amplitude(group_inds==grp,j)),1),std(squeeze(relative_MMN_amplitude(group_inds==grp,j)))/sqrt(sum(group_inds==grp)),'-s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1,'color','black')
            end
        else
            if amplitude_ranovatbl.pValueGG(find(strcmp(amplitude_ranovatbl.Row,'Diagnosis:Condition'))) < p_thresh && p<p_thresh %No longer Bonferroni as only do if diagnosis by region interaction is significant
                errorbar(grp-0.1,mean(squeeze(relative_MMN_amplitude(group_inds==grp,j)),1),std(squeeze(relative_MMN_amplitude(group_inds==grp,j)))/sqrt(sum(group_inds==grp)),'-s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red','LineWidth',1,'color','black')
            else
                errorbar(grp-0.1,mean(squeeze(relative_MMN_amplitude(group_inds==grp,j)),1),std(squeeze(relative_MMN_amplitude(group_inds==grp,j)))/sqrt(sum(group_inds==grp)),'-s','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor','black','LineWidth',1,'color','black')
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
saveas(relative_MMN_amplitude_plot,['./outputfigures/scalp/stats/Relative MMN amplitude.png']);


