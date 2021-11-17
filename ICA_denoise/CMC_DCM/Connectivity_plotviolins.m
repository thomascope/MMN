function Connectivity_plotviolins(dirname_DCM,diagnosis_list,source_names,thresh,Participant,p_data)
%A script for plotting the results of extDCM across all diagnoses by
%inter-regional connection
save_figures = 1;
PEB_focuses = {'A','B'};

all_folders = strsplit(dirname_DCM,'/');
all_folders = all_folders(~cellfun(@isempty,all_folders)); %For save location

Frequency_bands = [4 8; 8 20; 20 30; 30 45; 55 70]; % Avoid 50 because of electrical noise and filtering.
Frequence_band_names = {'Theta','Alpha','Beta','Low Gamma','High Gamma'};
post_thresh = 0.05; % Threshold for post-hoc tests;
num_perms = 5000; % Number of permutations for the by frequency tests
num_perms_t2 = 10; % Number of permutations for the t2 tests (not used, keep low as then this will be quick and never report as significant)

%Assign cell conditions as per spm_fx_cmc
cell_pops{1} = 'Superficial Pyramidal to Spiny Stellate';
cell_pops{2} = 'Superficial Pyramidal to Deep Pyramidal';
cell_pops{3} = 'Deep Pyramidal to Superficial Pyramidal';
cell_pops{4} = 'Deep Pyramidal to Inhibitory Interneuron';

addpath('/group/language/data/thomascope/MMN/ICA_denoise/Helperfiles')
thisdir = pwd;
mkdir('./CMC_figures/')

base_datapathstem = '/imaging/mlr/users/tc02/Holly_MMN/Coherence_Connectivity_Integrated_LOR/crosshem/'; %With sLORETA source reconstruction
addpath(['/group/language/data/thomascope/vespa/SPM12version/Standalone preprocessing pipeline/tc_source_stats/ojwoodford-export_fig-216b30e']);
addpath(['/group/language/data/thomascope/MMN/ICA_denoise/stdshade']);
addpath('/imaging/mlr/users/tc02/toolboxes/rsatoolbox/Engines/');

try
    set(0,'DefaultLegendAutoUpdate','off') %For more modern Matlab, which is too clever for its own good.
end

%groupstodo = {'Control' 'pca' 'bvFTD' 'nfvppa' 'MCI'};
for i = 1:length(Participant)
    all_diags{i} = Participant{i}.diag;
end

groupstodo = unique(all_diags,'stable');
cmap = colormap(parula(length(groupstodo)));

these_corr_sig_pairs = {1:2; 2:3; [2,4]; 3:4; 5:6; 6:7; [6,8]; 7:8; [1,5]; [2,6]; [3,7]; [4,8]};
sources = source_names;

averagesubtracted = 1;
highfreq = 1;
start_times = p_data.start_times;
end_times = p_data.end_times;
topfreqband = 40; % Include gamma
%topfreqband = 20;

analysis_type = {};
%analysis_type{end+1} = 'Granger'; %Influenced by SNR
%analysis_type{end+1} = 'icoh'; % Redundant given partialled versions
%analysis_type{end+1} = 'plv'; % Redundant given partialled versions
%analysis_type{end+1} = 'pdc'; % I don't trust this metric
analysis_type{end+1} = 'partial_plv';
analysis_type{end+1} = 'partial_icoh';

datapathstem = {};
filenames = {};
all_subjs = {};
all_granger_data = {};
all_mismatch_contrasts = {};
all_absolute_mismatch_contrasts = {};
all_random_granger_data = {};
all_random_mismatch_contrasts = {};
all_absolute_random_mismatch_contrasts = {};
demeaned_all_granger_data = {};
demeaned_all_mismatch_contrasts = {};

for i = 1:length(analysis_type)
    switch(analysis_type{i})
        case 'Granger'
            datapathstem{end+1} = [base_datapathstem 'granger/'];
        case 'icoh'
            datapathstem{end+1} = [base_datapathstem 'coh/'];
        case 'partial_icoh'
            datapathstem{end+1} = [base_datapathstem 'partial_coh/'];
        otherwise
            datapathstem{end+1} = [base_datapathstem analysis_type{i} '/'];
    end
    
    cd(datapathstem{end})
    
    filenames{end+1} = {};
    all_subjs{end+1} = [];
    
    runningtotal = 1;
    for j = 1:length(groupstodo)
        cd(groupstodo{j});
        thesefiles = dir(['*' num2str(start_times) '_' num2str(end_times) '_overall.mat']);
        all_subjs{i} = [all_subjs{i}; [repmat(j,size(thesefiles)), (runningtotal:(runningtotal+size(thesefiles,1)-1))']];
        runningtotal = runningtotal+size(thesefiles,1)';
        filenames{i} = [filenames{i}, thesefiles.name];
        cd(datapathstem{end})
    end
    
    group = all_subjs{i}(:,1)';
    
    s=1;
    load([datapathstem{i} groupstodo{group(s)} '/s' num2str(all_subjs{i}(s,2)) '_grangerdata_highfreq_averagesubtracted_100_' num2str(start_times) '_' num2str(end_times) '_overall'],'granger_data','foi');
    if topfreqband>length(foi)
        topfreqband = length(foi);
    end
    
    all_granger_data{i} = zeros(8,8,1,topfreqband,2,length(group));
    all_mismatch_contrasts{i} = zeros(8,8,1,topfreqband,length(group));
    all_absolute_mismatch_contrasts{i} = zeros(8,8,1,topfreqband,length(group));
    all_random_granger_data{i} = zeros(8,8,1,topfreqband,2,100,length(group));
    all_random_mismatch_contrasts{i} = zeros(8,8,1,topfreqband,100,length(group));
    all_absolute_random_mismatch_contrasts{i} = zeros(8,8,1,topfreqband,100,length(group));
    demeaned_all_granger_data{i} = zeros(8,8,1,topfreqband,2,length(group));
    demeaned_all_mismatch_contrasts{i} = zeros(8,8,1,topfreqband,length(group));
    
    for s = 1:length(group)
        load([datapathstem{i} groupstodo{group(s)} '/s' num2str(all_subjs{i}(s,2)) '_grangerdata_highfreq_averagesubtracted_100_' num2str(start_times) '_' num2str(end_times) '_overall'],'granger_data','foi');
        
        foi = foi(1:topfreqband);
        granger_data = granger_data(:,:,:,1:topfreqband,:,:);
        
        switch(analysis_type{i})
            case {'icoh', 'partial_icoh', 'plv', 'partial_plv'} % Not sure if taking the abs is correct for PLV theoretically, but I note that from/to reversals reverse the plv, so there must be a confound of directionality
                granger_data = abs(granger_data);
        end
        
        
        
        all_granger_data{i}(:,:,:,:,:,s) = granger_data(:,:,:,:,:,201);
        all_random_granger_data{i}(:,:,:,:,:,:,s) = granger_data(:,:,:,:,:,101:200); %Random permutations by trial
        demeaned_all_granger_data{i}(:,:,:,:,:,s) = granger_data(:,:,:,:,:,201) / mean(mean(mean(nanmean(granger_data(:,:,:,:,:,201)))));
        all_mismatch_contrasts{i}(:,:,:,:,s) = mean((granger_data(:,:,:,:,1,201)-granger_data(:,:,:,:,2,201)),5);
        all_absolute_mismatch_contrasts{i}(:,:,:,:,s) = mean((abs(granger_data(:,:,:,:,1,201))-abs(granger_data(:,:,:,:,2,201))),5);
        all_random_mismatch_contrasts{i}(:,:,:,:,:,s) = mean((granger_data(:,:,:,:,1,1:100)-granger_data(:,:,:,:,2,1:100)),5);
        all_absolute_random_mismatch_contrasts{i}(:,:,:,:,:,s) = mean((abs(granger_data(:,:,:,:,1,1:100))-abs(granger_data(:,:,:,:,2,1:100))),5);
        demeaned_all_mismatch_contrasts{i}(:,:,:,:,s) = mean((granger_data(:,:,:,:,1,201)-granger_data(:,:,:,:,2,201)),5) / mean(mean(mean(nanmean(granger_data(:,:,:,:,:,201)))));
    end
    
end

for this_focus = 1:length(PEB_focuses)
    PEB_focus = PEB_focuses{this_focus};
    cd([dirname_DCM 'PEB_secondlevel'])
    load(['PEB_' PEB_focus '_' cat(2,diagnosis_list{:}) '.mat'])
    cd(thisdir)
    template_PEB = PEB;
    
    if strcmp(PEB_focus,'A')  %Account for duplication of connections for NMDA and AMPA in A matrix
        BMA.Pp = BMA.Pp((strncmpi('A{1',BMA.Pnames,3)|strncmpi('A{3',BMA.Pnames,3)),:)+BMA.Pp((strncmpi('A{2',BMA.Pnames,3)|strncmpi('A{4',BMA.Pnames,3)),:);
        BMA.Pnames = BMA.Pnames(strncmpi('A{1',BMA.Pnames,3)|strncmpi('A{3',BMA.Pnames,3));
    end
    
    if strcmp(PEB_focus,'A')
        for this_contrast = 2:size(template_PEB.M.X,2)
            overall_main_effect_differences = zeros(size(template_PEB.M.X,2)-1,length(all_granger_data),length(BMA.Pnames),size(Frequency_bands,1));
            overall_null_main_effect_differences = zeros(size(template_PEB.M.X,2)-1,length(all_granger_data),length(BMA.Pnames),size(Frequency_bands,1),100);
            overall_null_main_effect_difference_means = zeros(size(template_PEB.M.X,2)-1,length(all_granger_data),length(BMA.Pnames),size(Frequency_bands,1));
            overall_null_main_effect_stdevs = zeros(size(template_PEB.M.X,2)-1,length(all_granger_data),length(BMA.Pnames),size(Frequency_bands,1));
            key_main_effect = zeros((length(BMA.Pnames)*length(all_granger_data)*(size(template_PEB.M.X,2)-1)),4);
        end
    elseif strcmp(PEB_focus,'B')
        for this_contrast = 2:size(template_PEB.M.X,2)
            overall_interaction_differences = zeros(size(template_PEB.M.X,2)-1,length(all_granger_data),length(BMA.Pnames),size(Frequency_bands,1));
            overall_null_interaction_differences = zeros(size(template_PEB.M.X,2)-1,length(all_granger_data),length(BMA.Pnames),size(Frequency_bands,1),100);
            overall_null_interaction_difference_means = zeros(size(template_PEB.M.X,2)-1,length(all_granger_data),length(BMA.Pnames),size(Frequency_bands,1));
            overall_null_interaction_stdevs = zeros(size(template_PEB.M.X,2)-1,length(all_granger_data),length(BMA.Pnames),size(Frequency_bands,1));
            key_interaction = zeros((length(BMA.Pnames)*length(all_granger_data)*(size(template_PEB.M.X,2)-1)),4);
        end
    end
    
    for this_loop = 1:(length(BMA.Pnames)*length(all_granger_data)*(size(template_PEB.M.X,2)-1))
        this_difference = ceil(this_loop/((size(template_PEB.M.X,2)-1)*length(all_granger_data)));
        this_measure = mod(ceil(this_loop/(size(template_PEB.M.X,2)-1))-1,length(all_granger_data))+1;
        this_contrast = 2+mod(this_loop-1,size(template_PEB.M.X,2)-1);
        disp(['Working on contrast ' num2str(this_contrast) ' measure ' num2str(this_measure) ' difference ' num2str(this_difference)])
        
        this_DCM_Pp_tf = BMA.Pp(this_difference,this_contrast)>thresh;
        
        if strcmp(PEB_focus,'A')
            key_main_effect(this_loop,:) = [this_contrast,this_measure,this_difference,this_DCM_Pp_tf];
        elseif strcmp(PEB_focus,'B')
            key_interaction(this_loop,:) = [this_contrast,this_measure,this_difference,this_DCM_Pp_tf];
        end
        this_connection = BMA.Pnames{this_difference};
        Connection_Split = strsplit(this_connection,[PEB_focus '{']);
        direction = Connection_Split{2}(1);
        if strcmp(direction,'1')
            direction = 'forwards';
        elseif strcmp(direction,'2')
            direction = 'backwards';
        end
        to = str2num(Connection_Split{2}(4));
        from = str2num(Connection_Split{2}(6));
        
        if to == from
            continue %No intrinsic connections in Granger/connectivity analysis.
        end
        
        if strcmp(PEB_focus,'A')
            for this_band = 1:size(Frequency_bands,1)
                these_fois = foi>Frequency_bands(this_band,1)&foi<Frequency_bands(this_band,2);
                overall_main_effect_differences(this_contrast,this_measure,this_difference,this_band) = mean(squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==1)),5)),1)))-mean(squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==-1)),5)),1)));
                overall_null_main_effect_differences(this_contrast,this_measure,this_difference,this_band,:) = mean(squeeze(mean(squeeze(mean(all_random_granger_data{this_measure}(from,to,:,these_fois,:,:,find(template_PEB.M.X(:,this_contrast)==1)),5)),1)),2)-mean(squeeze(mean(squeeze(mean(all_random_granger_data{this_measure}(from,to,:,these_fois,:,:,find(template_PEB.M.X(:,this_contrast)==-1)),5)),1)),2);
                overall_null_main_effect_stdevs(this_contrast,this_measure,this_difference,this_band) = std(mean(squeeze(mean(squeeze(mean(all_random_granger_data{this_measure}(from,to,:,these_fois,:,:,find(template_PEB.M.X(:,this_contrast)==1)),5)),1)),2)-mean(squeeze(mean(squeeze(mean(all_random_granger_data{this_measure}(from,to,:,these_fois,:,:,find(template_PEB.M.X(:,this_contrast)==-1)),5)),1)),2));
                overall_null_main_effect_difference_means(this_contrast,this_measure,this_difference,this_band) = mean(mean(squeeze(mean(squeeze(mean(all_random_granger_data{this_measure}(from,to,:,these_fois,:,:,find(template_PEB.M.X(:,this_contrast)==1)),5)),1)),2)-mean(squeeze(mean(squeeze(mean(all_random_granger_data{this_measure}(from,to,:,these_fois,:,:,find(template_PEB.M.X(:,this_contrast)==-1)),5)),1)),2));
            end
        elseif strcmp(PEB_focus,'B')
            for this_band = 1:size(Frequency_bands,1)
                these_fois = foi>Frequency_bands(this_band,1)&foi<Frequency_bands(this_band,2);
                overall_interaction_differences(this_contrast,this_measure,this_difference,this_band) = mean(squeeze(mean(squeeze(all_mismatch_contrasts{this_measure}(from,to,:,these_fois,find(template_PEB.M.X(:,this_contrast)==1))),1)))-mean(squeeze(mean(squeeze(all_mismatch_contrasts{this_measure}(from,to,:,these_fois,find(template_PEB.M.X(:,this_contrast)==-1))),1)));
                overall_null_interaction_differences(this_contrast,this_measure,this_difference,this_band,:) = mean(squeeze(mean(squeeze(all_random_mismatch_contrasts{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==1))),1)),2)-mean(squeeze(mean(squeeze(all_random_mismatch_contrasts{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==-1))),1)),2);
                overall_null_interaction_stdevs(this_contrast,this_measure,this_difference,this_band) = std(mean(squeeze(mean(squeeze(all_random_mismatch_contrasts{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==1))),1)),2)-mean(squeeze(mean(squeeze(all_random_mismatch_contrasts{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==-1))),1)),2));
                overall_null_interaction_difference_means(this_contrast,this_measure,this_difference,this_band) = mean(mean(squeeze(mean(squeeze(all_random_mismatch_contrasts{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==1))),1)),2)-mean(squeeze(mean(squeeze(all_random_mismatch_contrasts{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==-1))),1)),2));
            end
        end
    end
end

figure
set(gcf,'Position',[100 100 1600 800]);
set(gcf, 'PaperPositionMode', 'auto');
for contrast_of_interest = 2:5 % Controls against each patient group
    for this_measure = 1:length(all_granger_data)
        DCM_data_to_plot = [];
        non_DCM_data_to_plot = [];
        null_data_to_plot = [];
        null_data_stdevs = [];
        null_data_to_plot_all = [];
        all_data_to_plot = [];
        subplot(length(all_granger_data),4,(contrast_of_interest-1)+((this_measure-1)*4))
        for i = 1:size(key_main_effect,1)
            if ~(key_main_effect(i,1)==contrast_of_interest)
                continue
            elseif ~(key_main_effect(i,2)==this_measure)
                continue
            end
            if all(overall_main_effect_differences(key_main_effect(i,1),key_main_effect(i,2),key_main_effect(i,3),:)==0) %Intrinsic connection
                continue
            end
            disp(['Adding contrast ' num2str(key_main_effect(i,1)) ' measure ' num2str(key_main_effect(i,2)) ' difference ' num2str(key_main_effect(i,3))])
            if key_main_effect(i,4) == 1
                DCM_data_to_plot(end+1,:) = overall_main_effect_differences(key_main_effect(i,1),key_main_effect(i,2),key_main_effect(i,3),:);
            else
                non_DCM_data_to_plot(end+1,:) = overall_main_effect_differences(key_main_effect(i,1),key_main_effect(i,2),key_main_effect(i,3),:);
            end
            all_data_to_plot(end+1,:) = overall_main_effect_differences(key_main_effect(i,1),key_main_effect(i,2),key_main_effect(i,3),:);
            null_data_to_plot(end+1,:) = overall_null_main_effect_difference_means(key_main_effect(i,1),key_main_effect(i,2),key_main_effect(i,3),:);
            null_data_stdevs(end+1,:) = overall_null_main_effect_stdevs(key_main_effect(i,1),key_main_effect(i,2),key_main_effect(i,3),:);
            null_data_to_plot_all(end+1,:,:) = overall_null_main_effect_differences(key_main_effect(i,1),key_main_effect(i,2),key_main_effect(i,3),:,:);
        end
        hold on
        these_null_data = [];
        for this_band = 1:size(Frequency_bands,1)
            these_null_data(:,this_band) = reshape(null_data_to_plot_all(:,this_band,:),1,size(null_data_to_plot_all,1)*size(null_data_to_plot_all,3));
        end
        violin(these_null_data,'x',[1:5]-0.1)
        violin(all_data_to_plot,'x',[1:5]+0.1)
        
        xlim([0 6])
        set(gca,'xtick',[1:5],'xticklabels',Frequence_band_names,'XTickLabelRotation',45)
        data_bounds = [min(min(min(DCM_data_to_plot)),min(min(non_DCM_data_to_plot))),max(max(max(DCM_data_to_plot)),max(max(non_DCM_data_to_plot)))];
        ylim([data_bounds(1) - (diff(data_bounds)/10) data_bounds(2) + (diff(data_bounds)/6)])
        
        for this_band = 1:size(Frequency_bands,1)
            %[p(this_band), ~, ~] = permutationTest(these_null_data(:,this_band), all_data_to_plot(:,this_band), 2000);
            [p(this_band),h] = ranksum(all_data_to_plot(:,this_band),these_null_data(:,this_band));
            %[h p(this_band)] = ttest2(all_data_to_plot(:,this_band),these_null_data(:,this_band), 'vartype', 'unequal');
            if p(this_band) < 0.001
                text(this_band,data_bounds(2) + (diff(data_bounds)/10),'***')
            elseif p(this_band) < 0.01
                text(this_band,data_bounds(2) + (diff(data_bounds)/10),'**')
            elseif p(this_band) < 0.05
                text(this_band,data_bounds(2) + (diff(data_bounds)/10),'*')
            end
            %scatter((this_band)*ones(1,size(DCM_data_to_plot,1)),DCM_data_to_plot(:,this_band),12,'rx')
        end
        
        plot([0 6],[0,0],'k--')
        switch(analysis_type{this_measure})
            case 'Granger'
                title('Granger Causality')
            case 'icoh'
                title('Imaginary Coherence')
            case 'partial_icoh'
                title('Partial Imaginary Coherence')
            case 'plv'
                title('Phase Locking Value')
            case 'partial_plv'
                title('Partial Phase Locking Value')
            otherwise
                title(analysis_type{this_measure},'interpreter','none');
        end
    end
    %     suptitle([Participant{min(find(template_PEB.M.X(:,contrast_of_interest)==1))}.diag ' minus ' Participant{min(find(template_PEB.M.X(:,contrast_of_interest)==-1))}.diag ' Main effect'])
end
savestring = ['./CMC_figures/Connectivity_Violin_' num2str(start_times) '_' num2str(end_times) '.pdf'];
if save_figures
    drawnow
    saveas(gcf,savestring)
    eval(['export_fig ' savestring(1:end-3) 'png']);
end

figure
set(gcf,'Position',[100 100 1600 800]);
set(gcf, 'PaperPositionMode', 'auto');
for contrast_of_interest = 2:5 % Controls against each patient group
    for this_measure = 1:length(all_granger_data)
        DCM_data_to_plot = [];
        non_DCM_data_to_plot = [];
        null_data_to_plot = [];
        null_data_stdevs = [];
        null_data_to_plot_all = [];
        all_data_to_plot = [];
        subplot(length(all_granger_data),4,(contrast_of_interest-1)+((this_measure-1)*4))
        for i = 1:size(key_main_effect,1)
            if ~(key_main_effect(i,1)==contrast_of_interest)
                continue
            elseif ~(key_main_effect(i,2)==this_measure)
                continue
            end
            if all(overall_main_effect_differences(key_main_effect(i,1),key_main_effect(i,2),key_main_effect(i,3),:)==0) %Intrinsic connection
                continue
            end
            disp(['Adding contrast ' num2str(key_main_effect(i,1)) ' measure ' num2str(key_main_effect(i,2)) ' difference ' num2str(key_main_effect(i,3))])
            if key_main_effect(i,4) == 1
                DCM_data_to_plot(end+1,:) = overall_main_effect_differences(key_main_effect(i,1),key_main_effect(i,2),key_main_effect(i,3),:);
            else
                non_DCM_data_to_plot(end+1,:) = overall_main_effect_differences(key_main_effect(i,1),key_main_effect(i,2),key_main_effect(i,3),:);
            end
            all_data_to_plot(end+1,:) = overall_main_effect_differences(key_main_effect(i,1),key_main_effect(i,2),key_main_effect(i,3),:);
            null_data_to_plot(end+1,:) = overall_null_main_effect_difference_means(key_main_effect(i,1),key_main_effect(i,2),key_main_effect(i,3),:);
            null_data_stdevs(end+1,:) = overall_null_main_effect_stdevs(key_main_effect(i,1),key_main_effect(i,2),key_main_effect(i,3),:);
            null_data_to_plot_all(end+1,:,:) = overall_null_main_effect_differences(key_main_effect(i,1),key_main_effect(i,2),key_main_effect(i,3),:,:);
        end
        hold on
        these_null_data = [];
        for this_band = 1:size(Frequency_bands,1)
            these_null_data(:,this_band) = reshape(null_data_to_plot_all(:,this_band,:),1,size(null_data_to_plot_all,1)*size(null_data_to_plot_all,3));
        end
        violin(these_null_data,'x',[1:5]-0.1)
        violin(all_data_to_plot,'x',[1:5]+0.1)
        
        xlim([0 6])
        set(gca,'xtick',[1:5],'xticklabels',Frequence_band_names,'XTickLabelRotation',45)
        data_bounds = [min(min(min(DCM_data_to_plot)),min(min(non_DCM_data_to_plot))),max(max(max(DCM_data_to_plot)),max(max(non_DCM_data_to_plot)))];
        ylim([data_bounds(1) - (diff(data_bounds)/10) data_bounds(2) + (diff(data_bounds)/6)])
        
        for this_band = 1:size(Frequency_bands,1)
            %[p(this_band), ~, ~] = permutationTest(these_null_data(:,this_band), all_data_to_plot(:,this_band), 2000);
            %[h p(this_band)] = ttest2(all_data_to_plot(:,this_band),these_null_data(:,this_band), 'vartype', 'unequal');
            [p(this_band),h] = ranksum(all_data_to_plot(:,this_band),these_null_data(:,this_band));
            p(this_band) = p(this_band)*5; %Bonferroni correct across frequency bands
            if p(this_band) < 0.001
                text(this_band,data_bounds(2) + (diff(data_bounds)/10),'***')
            elseif p(this_band) < 0.01
                text(this_band,data_bounds(2) + (diff(data_bounds)/10),'**')
            elseif p(this_band) < 0.05
                text(this_band,data_bounds(2) + (diff(data_bounds)/10),'*')
            end
            %scatter((this_band)*ones(1,size(DCM_data_to_plot,1)),DCM_data_to_plot(:,this_band),12,'rx')
        end
        
        plot([0 6],[0,0],'k--')
        switch(analysis_type{this_measure})
            case 'Granger'
                title('Granger Causality')
            case 'icoh'
                title('Imaginary Coherence')
            case 'partial_icoh'
                title('Partial Imaginary Coherence')
            case 'plv'
                title('Phase Locking Value')
            case 'partial_plv'
                title('Partial Phase Locking Value')
            otherwise
                title(analysis_type{this_measure},'interpreter','none');
        end
    end
    %     suptitle([Participant{min(find(template_PEB.M.X(:,contrast_of_interest)==1))}.diag ' minus ' Participant{min(find(template_PEB.M.X(:,contrast_of_interest)==-1))}.diag ' Main effect'])
end
savestring = ['./CMC_figures/Connectivity_Violin_Bonferroni_' num2str(start_times) '_' num2str(end_times) '.pdf'];
if save_figures
    drawnow
    saveas(gcf,savestring)
    eval(['export_fig ' savestring(1:end-3) 'png']);
end

%Now break down by connection type
Auditory = {'2,1';'1,2';'5,6';'6,5'};
Cross_Hem = {'7,3';'3,7';'8,4';'4,8';'2,6';'6,2'};
Aud_MD = {'3,2';'2,3';'2,4';'4,2';'7,6';'6,7';'6,8';'8,6'};
MD = {'3,4';'4,3';'7,8';'8,7'};

this_focus = 1;
PEB_focus = PEB_focuses{this_focus};
cd([dirname_DCM 'PEB_secondlevel'])
load(['PEB_' PEB_focus '_' cat(2,diagnosis_list{:}) '.mat'])
cd(thisdir)
template_PEB = PEB;

if strcmp(PEB_focus,'A')  %Account for duplication of connections for NMDA and AMPA in A matrix
    BMA.Pp = BMA.Pp((strncmpi('A{1',BMA.Pnames,3)|strncmpi('A{3',BMA.Pnames,3)),:)+BMA.Pp((strncmpi('A{2',BMA.Pnames,3)|strncmpi('A{4',BMA.Pnames,3)),:);
    BMA.Pnames = BMA.Pnames(strncmpi('A{1',BMA.Pnames,3)|strncmpi('A{3',BMA.Pnames,3));
end

figure
set(gcf,'Position',[100 100 1600 800]);
set(gcf, 'PaperPositionMode', 'auto');
all_null_data = [];
all_Auditory_data = [];
all_Cross_Hem_data = [];
all_Aud_MD_data = [];
all_MD_data = [];

for contrast_of_interest = 2:5 % Controls against each patient group
    for this_measure = 1:length(all_granger_data)
        Auditory_data_to_plot = [];
        Cross_Hem_data_to_plot = [];
        Aud_MD_data_to_plot = [];
        MD_data_to_plot = [];
        null_data_to_plot = [];
        null_data_stdevs = [];
        null_data_to_plot_all = [];
        all_data_to_plot = [];
        subplot(length(all_granger_data),4,(contrast_of_interest-1)+((this_measure-1)*4))
        
        for i = 1:size(key_main_effect,1)
            if ~(key_main_effect(i,1)==contrast_of_interest)
                continue
            elseif ~(key_main_effect(i,2)==this_measure)
                continue
            end
            if all(overall_main_effect_differences(key_main_effect(i,1),key_main_effect(i,2),key_main_effect(i,3),:)==0) %Intrinsic connection
                continue
            end
            disp(['Adding contrast ' num2str(key_main_effect(i,1)) ' measure ' num2str(key_main_effect(i,2)) ' difference ' num2str(key_main_effect(i,3))])
            
            this_connection = BMA.Pnames{key_main_effect(i,3)};
            Connection_Split = strsplit(this_connection,[PEB_focus '{']);
            this_connection = Connection_Split{2}(4:6);
            
            if any(strcmp(this_connection,Auditory))
                Auditory_data_to_plot(end+1,:) = overall_main_effect_differences(key_main_effect(i,1),key_main_effect(i,2),key_main_effect(i,3),:);
            elseif any(strcmp(this_connection,Cross_Hem))
                Cross_Hem_data_to_plot(end+1,:) = overall_main_effect_differences(key_main_effect(i,1),key_main_effect(i,2),key_main_effect(i,3),:);
            elseif any(strcmp(this_connection,Aud_MD))
                Aud_MD_data_to_plot(end+1,:) = overall_main_effect_differences(key_main_effect(i,1),key_main_effect(i,2),key_main_effect(i,3),:);
            elseif any(strcmp(this_connection,MD))
                MD_data_to_plot(end+1,:) = overall_main_effect_differences(key_main_effect(i,1),key_main_effect(i,2),key_main_effect(i,3),:);
            else
                error([this_connection ' not found'])
            end
            all_data_to_plot(end+1,:) = overall_main_effect_differences(key_main_effect(i,1),key_main_effect(i,2),key_main_effect(i,3),:);
            null_data_to_plot(end+1,:) = overall_null_main_effect_difference_means(key_main_effect(i,1),key_main_effect(i,2),key_main_effect(i,3),:);
            null_data_stdevs(end+1,:) = overall_null_main_effect_stdevs(key_main_effect(i,1),key_main_effect(i,2),key_main_effect(i,3),:);
            null_data_to_plot_all(end+1,:,:) = overall_null_main_effect_differences(key_main_effect(i,1),key_main_effect(i,2),key_main_effect(i,3),:,:);
        end
        hold on
        these_null_data = [];
        for this_band = 1:size(Frequency_bands,1)
            these_null_data(:,this_band) = reshape(null_data_to_plot_all(:,this_band,:),1,size(null_data_to_plot_all,1)*size(null_data_to_plot_all,3));
        end
        violin(these_null_data,'x',[1:5]-0.3)
        violin(Auditory_data_to_plot,'x',[1:5]-0.15)
        violin(Cross_Hem_data_to_plot,'x',[1:5])
        violin(Aud_MD_data_to_plot,'x',[1:5]+0.15)
        violin(MD_data_to_plot,'x',[1:5]+0.3)
        
        xlim([0 6])
        set(gca,'xtick',[1:5],'xticklabels',Frequence_band_names,'XTickLabelRotation',45)
        data_bounds = [min(min(all_data_to_plot)),max(max(all_data_to_plot))];
        ylim([data_bounds(1) - (diff(data_bounds)/10) data_bounds(2) + (diff(data_bounds)/6)])
        
        %         for this_band = 1:size(Frequency_bands,1)
        %             %[p(this_band), ~, ~] = permutationTest(these_null_data(:,this_band), all_data_to_plot(:,this_band), 2000);
        %             [h p(this_band)] = ttest2(all_data_to_plot(:,this_band),these_null_data(:,this_band), 'vartype', 'unequal');
        %             if p(this_band) < 0.001
        %                 text(this_band,data_bounds(2) + (diff(data_bounds)/10),'***')
        %             elseif p(this_band) < 0.01
        %                 text(this_band,data_bounds(2) + (diff(data_bounds)/10),'**')
        %             elseif p(this_band) < 0.05
        %                 text(this_band,data_bounds(2) + (diff(data_bounds)/10),'*')
        %             end
        %             scatter((this_band)*ones(1,size(DCM_data_to_plot,1)),DCM_data_to_plot(:,this_band),12,'rx')
        %         end
        
        plot([0 6],[0,0],'k--')
        switch(analysis_type{this_measure})
            case 'Granger'
                title('Granger Causality')
            case 'icoh'
                title('Imaginary Coherence')
            case 'partial_icoh'
                title('Partial Imaginary Coherence')
            case 'plv'
                title('Phase Locking Value')
            case 'partial_plv'
                title('Partial Phase Locking Value')
            otherwise
                title(analysis_type{this_measure},'interpreter','none');
        end
        all_null_data(contrast_of_interest,this_measure,:,:) = these_null_data;
        all_Auditory_data(contrast_of_interest,this_measure,:,:) = Auditory_data_to_plot;
        all_Cross_Hem_data(contrast_of_interest,this_measure,:,:) = Cross_Hem_data_to_plot;
        all_Aud_MD_data(contrast_of_interest,this_measure,:,:) = Aud_MD_data_to_plot;
        all_MD_data(contrast_of_interest,this_measure,:,:) = MD_data_to_plot;
    end
    %     suptitle([Participant{min(find(template_PEB.M.X(:,contrast_of_interest)==1))}.diag ' minus ' Participant{min(find(template_PEB.M.X(:,contrast_of_interest)==-1))}.diag ' Main effect']
end

for this_measure = 1:length(all_granger_data)
    figure
    set(gcf,'Position',[100 100 1600 800]);
    set(gcf, 'PaperPositionMode', 'auto');
    scatter_colors = [0,0,0;1,0,0];
    for this_band = 1:size(Frequency_bands,1)
        subplot(4,size(Frequency_bands,1),this_band+((1-1)*size(Frequency_bands,1)))
        violin(reshape((squeeze(all_null_data(2:end,this_measure,:,this_band))'),numel((squeeze(all_null_data(2:end,this_measure,:,this_band))')),1))
        hold on
        for this_group = 2:5
            these_colors = scatter_colors((abs(unique(squeeze(all_Auditory_data(this_group,this_measure,:,this_band))))>max(abs(reshape((squeeze(all_null_data(2:end,this_measure,:,this_band))'),numel((squeeze(all_null_data(2:end,this_measure,:,this_band))')),1))))+1,:);
            scatter(repmat(this_group,size(unique(squeeze(all_Auditory_data(this_group,this_measure,:,this_band))),1),1),unique(squeeze(all_Auditory_data(this_group,this_measure,:,this_band))),16,these_colors,'filled')
        end
        data_bounds = [min(min(squeeze(all_Auditory_data(2:end,this_measure,:,this_band))')),max(max(squeeze(all_Auditory_data(2:end,this_measure,:,this_band))'))];
        ylim([data_bounds(1) - (diff(data_bounds)/10) data_bounds(2) + (diff(data_bounds)/6)])
        %violin(squeeze(all_Auditory_data(2:end,this_measure,:,this_band))','x',[2:5])
        set(gca,'xtick',[1:5],'xticklabels',{'null','bvFTD','nfvPPA','PCA','ADMCI'},'XTickLabelRotation',45)
        xlim([0 6])
        plot([0 6],[0,0],'k--')
        title(['Auditory ' Frequence_band_names{this_band}])
        
        subplot(4,size(Frequency_bands,1),this_band+((2-1)*size(Frequency_bands,1)))
        violin(reshape((squeeze(all_null_data(2:end,this_measure,:,this_band))'),numel((squeeze(all_null_data(2:end,this_measure,:,this_band))')),1))
        hold on
        for this_group = 2:5
            these_colors = scatter_colors((abs(unique(squeeze(all_Cross_Hem_data(this_group,this_measure,:,this_band))))>max(abs(reshape((squeeze(all_null_data(2:end,this_measure,:,this_band))'),numel((squeeze(all_null_data(2:end,this_measure,:,this_band))')),1))))+1,:);
            scatter(repmat(this_group,size(unique(squeeze(all_Cross_Hem_data(this_group,this_measure,:,this_band))),1),1),unique(squeeze(all_Cross_Hem_data(this_group,this_measure,:,this_band))),16,these_colors,'filled')
        end
        data_bounds = [min(min(squeeze(all_Cross_Hem_data(2:end,this_measure,:,this_band))')),max(max(squeeze(all_Cross_Hem_data(2:end,this_measure,:,this_band))'))];
        ylim([data_bounds(1) - (diff(data_bounds)/10) data_bounds(2) + (diff(data_bounds)/6)])
        %violin(squeeze(all_Cross_Hem_data(2:end,this_measure,:,this_band))','x',[2:5])
        set(gca,'xtick',[1:5],'xticklabels',{'null','bvFTD','nfvPPA','PCA','ADMCI'},'XTickLabelRotation',45)
        xlim([0 6])
        plot([0 6],[0,0],'k--')
        title(['Cross-Hem ' Frequence_band_names{this_band}])
        
        subplot(4,size(Frequency_bands,1),this_band+((3-1)*size(Frequency_bands,1)))
        violin(reshape((squeeze(all_null_data(2:end,this_measure,:,this_band))'),numel((squeeze(all_null_data(2:end,this_measure,:,this_band))')),1))
        hold on
        for this_group = 2:5
            these_colors = scatter_colors((abs(unique(squeeze(all_Aud_MD_data(this_group,this_measure,:,this_band))))>max(abs(reshape((squeeze(all_null_data(2:end,this_measure,:,this_band))'),numel((squeeze(all_null_data(2:end,this_measure,:,this_band))')),1))))+1,:);
            scatter(repmat(this_group,size(unique(squeeze(all_Aud_MD_data(this_group,this_measure,:,this_band))),1),1),unique(squeeze(all_Aud_MD_data(this_group,this_measure,:,this_band))),16,these_colors,'filled')
        end
        data_bounds = [min(min(squeeze(all_Aud_MD_data(2:end,this_measure,:,this_band))')),max(max(squeeze(all_Aud_MD_data(2:end,this_measure,:,this_band))'))];
        ylim([data_bounds(1) - (diff(data_bounds)/10) data_bounds(2) + (diff(data_bounds)/6)])
        %violin(squeeze(all_Aud_MD_data(2:end,this_measure,:,this_band))','x',[2:5])
        set(gca,'xtick',[1:5],'xticklabels',{'null','bvFTD','nfvPPA','PCA','ADMCI'},'XTickLabelRotation',45)
        xlim([0 6])
        plot([0 6],[0,0],'k--')
        title(['Aud-MD ' Frequence_band_names{this_band}])
        
        subplot(4,size(Frequency_bands,1),this_band+((4-1)*size(Frequency_bands,1)))
        violin(reshape((squeeze(all_null_data(2:end,this_measure,:,this_band))'),numel((squeeze(all_null_data(2:end,this_measure,:,this_band))')),1))
        hold on
        for this_group = 2:5
            these_colors = scatter_colors((abs(unique(squeeze(all_MD_data(this_group,this_measure,:,this_band))))>max(abs(reshape((squeeze(all_null_data(2:end,this_measure,:,this_band))'),numel((squeeze(all_null_data(2:end,this_measure,:,this_band))')),1))))+1,:);
            scatter(repmat(this_group,size(unique(squeeze(all_MD_data(this_group,this_measure,:,this_band))),1),1),unique(squeeze(all_MD_data(this_group,this_measure,:,this_band))),16,these_colors,'filled')
        end
        data_bounds = [min(min(squeeze(all_MD_data(2:end,this_measure,:,this_band))')),max(max(squeeze(all_MD_data(2:end,this_measure,:,this_band))'))];
        ylim([data_bounds(1) - (diff(data_bounds)/10) data_bounds(2) + (diff(data_bounds)/6)])
        %violin(squeeze(all_MD_data(2:end,this_measure,:,this_band))','x',[2:5])
        set(gca,'xtick',[1:5],'xticklabels',{'null','bvFTD','nfvPPA','PCA','ADMCI'},'XTickLabelRotation',45)
        xlim([0 6])
        plot([0 6],[0,0],'k--')
        title(['Multiple Demand ' Frequence_band_names{this_band}])
    end
    
    
    switch(analysis_type{this_measure})
        case 'Granger'
            suptitle('Granger Causality')
        case 'icoh'
            suptitle('Imaginary Coherence')
        case 'partial_icoh'
            suptitle('Partial Imaginary Coherence')
        case 'plv'
            suptitle('Phase Locking Value')
        case 'partial_plv'
            suptitle('Partial Phase Locking Value')
        otherwise
            suptitle(analysis_type{this_measure},'interpreter','none');
    end
    savestring = ['./CMC_figures/Connectivity_Byconnection_' analysis_type{this_measure} '_' num2str(start_times) '_' num2str(end_times) '.pdf'];
    if save_figures
        drawnow
        saveas(gcf,savestring)
        eval(['export_fig ' savestring(1:end-3) 'png']);
    end
end

