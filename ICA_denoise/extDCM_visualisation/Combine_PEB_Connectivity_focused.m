function Combine_PEB_Connectivity_focused(dirname_DCM,diagnosis_list,source_names,thresh,Participant)
%A script for plotting the results of extDCM across all diagnoses by
%inter-regional connection
Frequency_bands = [4 8; 8 20; 20 30; 30 45; 55 70]; % Avoid 50 because of electrical noise and filtering.
Frequence_band_names = {'Theta','Alpha','Beta','Low Gamma','High Gamma'};

addpath('/group/language/data/thomascope/MMN/ICA_denoise/Helperfiles')
thisdir = pwd;
mkdir([thisdir '/circuit_diagrams'])
cd([dirname_DCM 'PEB_secondlevel'])

load('PEB_A_Overall.mat')
template_PEB = PEB_Overall;

assert(all(template_PEB.M.X(:,1)==1),'The first column of the PEB of PEBs contrast should be all ones, check please.')

base_datapathstem = '/imaging/tc02/Holly_MMN/Coherence_Connectivity_Integrated_LOR/crosshem/'; %With sLORETA source reconstruction
addpath(['/group/language/data/thomascope/vespa/SPM12version/Standalone preprocessing pipeline/tc_source_stats/ojwoodford-export_fig-216b30e']);
addpath(['/group/language/data/thomascope/MMN/ICA_denoise/stdshade']);
addpath('/imaging/tc02/toolboxes/rsatoolbox/Engines/');

set(0,'DefaultLegendAutoUpdate','off')

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
start_times = 0;
end_times = 500;
topfreqband = 40; % Include gamma
%topfreqband = 20;

save_figures = 1;

analysis_type = {};
analysis_type{end+1} = 'Granger';
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
        thesefiles = dir('*overall*mat');
        all_subjs{i} = [all_subjs{i}; [repmat(j,size(thesefiles)), (runningtotal:(runningtotal+size(thesefiles,1)-1))']];
        runningtotal = runningtotal+size(thesefiles,1)';
        filenames{i} = [filenames{i}, thesefiles.name];
        cd(datapathstem{end})
    end
    
    group = all_subjs{i}(:,1)';
    
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

cd(thisdir)

for this_contrast = 2:size(template_PEB.M.X,2)
    these_differences = find(BMA_Overall.Pp(:,this_contrast)>thresh);
    
    for this_difference = 1:length(these_differences)
        
        this_connection = BMA_Overall.Pnames{these_differences(this_difference)};
        Condition_Split = strsplit(this_connection,'Covariate ');
        condition = str2num(Condition_Split{2}(1));
        
        Connection_Split = strsplit(this_connection,'A{');
        direction = Connection_Split{2}(1);
        if strcmp(direction,'1')
            direction = 'forwards';
        elseif strcmp(direction,'2')
            direction = 'backwards';
        end
        to = str2num(Connection_Split{2}(4));
        from = str2num(Connection_Split{2}(6));
        
        
        if condition == 1
            if BMA_Overall.Ep(these_differences(this_difference),this_contrast)>0
                disp([diagnosis_list{find(template_PEB.M.X(:,this_contrast)==1)} ' stronger than ' diagnosis_list{find(template_PEB.M.X(:,this_contrast)==-1)} ' for ' direction ' ' source_names{from} ' to ' source_names{to} ' with posterior probability ' num2str(BMA_Overall.Pp(these_differences(this_difference),this_contrast))])
            else
                disp([diagnosis_list{find(template_PEB.M.X(:,this_contrast)==-1)} ' stronger than ' diagnosis_list{find(template_PEB.M.X(:,this_contrast)==1)} ' for ' direction ' ' source_names{from} ' to ' source_names{to} ' with posterior probability ' num2str(BMA_Overall.Pp(these_differences(this_difference),this_contrast))])
            end
        else
            disp(['Interaction between ' diagnosis_list{find(template_PEB.M.X(:,this_contrast)==-1)} ' and ' diagnosis_list{find(template_PEB.M.X(:,this_contrast)==1)} ' for ' direction ' ' source_names{from} ' to ' source_names{to} ' with posterior probability ' num2str(BMA_Overall.Pp(these_differences(this_difference),this_contrast))])
        end
        
        this_fig = figure;
        set(gcf,'Position',[100 100 1600 800]);
        set(gcf, 'PaperPositionMode', 'auto');
        clear linehandle legend_text
        % Reminder: size(all_granger_data{1}) =     8     8     1    20     2   123
        % From, to, time_window, foi, STD-DVT, subj
        for this_subplot = 1:4*length(all_granger_data)
            subplot(4,length(all_granger_data),this_subplot)
            this_measure = mod(this_subplot,length(all_granger_data));
            if this_measure == 0;  this_measure = length(all_granger_data); end
            this_connectivity_contrast = ceil(this_subplot/length(all_granger_data));
            
            if this_connectivity_contrast == 1
                % First look at the overall directionality by group - de-meaned
                % by individual, not meaningful for the icoh contrasts
                difference_demeaned = demeaned_all_granger_data{this_measure}(from,to,:,:,:,:)-demeaned_all_granger_data{this_measure}(to,from,:,:,:,:);
                hold on
                plot([0 ceil(max(foi))],[0 0],'k--','LineWidth',1);
                xlim([0 ceil(max(foi)/10)*10])
                legend_text = cell(1,2*length(groupstodo));
                for grp = 1:2:2*length(groupstodo)
                    if ceil(grp/2) == find(template_PEB.M.X(:,this_contrast)==1) || ceil(grp/2) == find(template_PEB.M.X(:,this_contrast)==-1)
                        linehandle(grp) = stdshade_TEC_cmap(squeeze(mean(difference_demeaned(1,1,:,:,:,group==ceil(grp/2)),5))',0.2,cmap(ceil(grp/2),:),foi,1,1,':');
                        legend_text{grp} = [groupstodo{ceil(grp/2)}];
                    end
                end
                if this_measure == length(all_granger_data)
                    legend(linehandle(~cellfun('isempty',legend_text)),legend_text(~cellfun('isempty',legend_text)));
                end
                title(['Directionality ' analysis_type{this_measure}], 'Interpreter', 'none')
                xlabel('Frequency, Hz')
            elseif this_connectivity_contrast == 2
                % Next look at the overall connection strength by group -
                % confounded by SNR for Granger and perhaps pdc
                hold on
                plot([0 ceil(max(foi))],[0 0],'k--','LineWidth',1);
                xlim([0 ceil(max(foi)/10)*10])
                for grp = 1:2:2*length(groupstodo)
                    if ceil(grp/2) == find(template_PEB.M.X(:,this_contrast)==1) || ceil(grp/2) == find(template_PEB.M.X(:,this_contrast)==-1)
                        stdshade_TEC_cmap(squeeze(mean(all_granger_data{this_measure}(from,to,:,:,:,group==ceil(grp/2)),5))',0.2,cmap(ceil(grp/2),:),foi,1,1,':');
                    end
                end
                if condition == 1
                    for this_band = 1:size(Frequency_bands,1)
                        these_fois = foi>Frequency_bands(this_band,1)&foi<Frequency_bands(this_band,2);
                        [~, pval] = ttest2(squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,group==find(template_PEB.M.X(:,this_contrast)==1)),5)),1)),squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,group==find(template_PEB.M.X(:,this_contrast)==-1)),5)),1)),'vartype','unequal');
                        if pval <0.05
                            if mean(squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,group==find(template_PEB.M.X(:,this_contrast)==1)),5)),1)))>mean(squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,group==find(template_PEB.M.X(:,this_contrast)==-1)),5)),1)))
                                disp([groupstodo{find(template_PEB.M.X(:,this_contrast)==1)} ' stronger than ' groupstodo{find(template_PEB.M.X(:,this_contrast)==-1)} ' in the ' Frequence_band_names{this_band} ' band using metric ' analysis_type{this_measure} ' p= ' num2str(pval)])
                            else
                                disp([groupstodo{find(template_PEB.M.X(:,this_contrast)==-1)} ' stronger than ' groupstodo{find(template_PEB.M.X(:,this_contrast)==1)} ' in the ' Frequence_band_names{this_band} ' band using metric ' analysis_type{this_measure} ' p= ' num2str(pval)])
                            end
                        end
                    end
                end
                title(['Strength ' analysis_type{this_measure}], 'Interpreter', 'none')
                xlabel('Frequency, Hz')
            elseif this_connectivity_contrast == 3
                % Next look at the modulation by STD-DVT
                hold on
                plot([0 ceil(max(foi))],[0 0],'k--','LineWidth',1);
                xlim([0 ceil(max(foi)/10)*10])
                for grp = 1:2:2*length(groupstodo)
                    if ceil(grp/2) == find(template_PEB.M.X(:,this_contrast)==1) || ceil(grp/2) == find(template_PEB.M.X(:,this_contrast)==-1)
                        stdshade_TEC_cmap(squeeze(all_mismatch_contrasts{this_measure}(from,to,:,:,group==ceil(grp/2)))',0.2,cmap(ceil(grp/2),:),foi,1,1,':');
                    end
                end
                if condition ~= 1
                    for this_band = 1:size(Frequency_bands,1)
                        these_fois = foi>Frequency_bands(this_band,1)&foi<Frequency_bands(this_band,2);
                        [~, pval] = ttest2(squeeze(mean(squeeze(all_mismatch_contrasts{this_measure}(from,to,:,these_fois,group==find(template_PEB.M.X(:,this_contrast)==1))),1)),squeeze(mean(squeeze(all_mismatch_contrasts{this_measure}(from,to,:,these_fois,group==find(template_PEB.M.X(:,this_contrast)==-1))),1)),'vartype','unequal');
                        if pval <0.05
                            disp([groupstodo{find(template_PEB.M.X(:,this_contrast)==1)} ' interacts with ' groupstodo{find(template_PEB.M.X(:,this_contrast)==-1)} ' in the ' Frequence_band_names{this_band} ' band using metric ' analysis_type{this_measure} ' p= ' num2str(pval)])
                        end
                    end
                end
                title(['STD-DVT ' analysis_type{this_measure}], 'Interpreter', 'none')
                xlabel('Frequency, Hz')
            elseif this_connectivity_contrast == 4
                % Finally look at the power for the reverse direction of the same pair -
                % only meaningful for some measures
                % Next look at the overall connection strength by group -
                % confounded by SNR for Granger and perhaps pdc
                hold on
                plot([0 ceil(max(foi))],[0 0],'k--','LineWidth',1);
                xlim([0 ceil(max(foi)/10)*10])
                for grp = 1:2:2*length(groupstodo)
                    if ceil(grp/2) == find(template_PEB.M.X(:,this_contrast)==1) || ceil(grp/2) == find(template_PEB.M.X(:,this_contrast)==-1)
                        stdshade_TEC_cmap(squeeze(mean(all_granger_data{this_measure}(to,from,:,:,:,group==ceil(grp/2)),5))',0.2,cmap(ceil(grp/2),:),foi,1,1,':');
                    end
                end
                title(['Strength ' analysis_type{this_measure}], 'Interpreter', 'none')
                xlabel('Frequency, Hz')
                
            end
        end
        suptitle(['Connectivity from ' source_names{from}  ' to ' source_names{to}]);
        if ceil(max(foi)) > 45
            savestring = ['./figures/' source_names{from} '_' source_names{to} '_Multifig_highfreq_focused.pdf'];
        else
            savestring = ['./figures/' source_names{from} '_' source_names{to} '_Multifig_focused.pdf'];
        end
        
        savestring = strrep(savestring,' ','_');
        if save_figures
            print(savestring,'-depsc','-painters'); %eval(['export_fig ' savestring ' -transparent']);
            eval(['export_fig ' savestring(1:end-3) 'png -transparent']);
        end
        
    end
end

cd([dirname_DCM 'PEB_secondlevel'])
load('PEB_D_Overall.mat')
cd(thisdir)

for this_contrast = 2:size(template_PEB.M.X,2)
    these_differences = find(BMA_Overall.Pp(:,this_contrast)>thresh);
    
    for this_difference = 1:length(these_differences)
        
        delays = {'local','cortico-cortical','cor-thalamo-cor'};
        
        this_connection = BMA_Overall.Pnames{these_differences(this_difference)};
        Condition_Split = strsplit(this_connection,'Covariate ');
        condition = str2num(Condition_Split{2}(1));
        
        Connection_Split = strsplit(this_connection,'D(');
        connection = delays{str2num(Connection_Split{2}(1))};
        
        if condition == 1
            if BMA_Overall.Ep(these_differences(this_difference),this_contrast)>0
                disp([diagnosis_list{find(template_PEB.M.X(:,this_contrast)==1)} ' longer than ' diagnosis_list{find(template_PEB.M.X(:,this_contrast)==-1)} ' for ' connection ' with posterior probability ' num2str(BMA_Overall.Pp(these_differences(this_difference),this_contrast))])
            else
                disp([diagnosis_list{find(template_PEB.M.X(:,this_contrast)==-1)} ' longer than ' diagnosis_list{find(template_PEB.M.X(:,this_contrast)==1)} ' for ' connection ' with posterior probability ' num2str(BMA_Overall.Pp(these_differences(this_difference),this_contrast))])
            end
        else
            disp(['Interaction between ' diagnosis_list{find(template_PEB.M.X(:,this_contrast)==-1)} ' and ' diagnosis_list{find(template_PEB.M.X(:,this_contrast)==1)} ' for ' connection ' with posterior probability ' num2str(BMA_Overall.Pp(these_differences(this_difference),this_contrast))])
        end
        
    end
end


cd([dirname_DCM 'PEB_secondlevel'])
load('PEB_A_Overall_combined.mat')
template_PEB = PEB_Overall;
cd(thisdir)
new_groupstodo = {'Control','All_FTD','All_AD'};

assert(all(template_PEB.M.X(:,1)==1),'The first column of the PEB of PEBs contrast should be all ones, check please.')

for this_contrast = 2:size(template_PEB.M.X,2)
    these_differences = find(BMA_Overall.Pp(:,this_contrast)>thresh);
    
    for this_difference = 1:length(these_differences)
        
        this_connection = BMA_Overall.Pnames{these_differences(this_difference)};
        Condition_Split = strsplit(this_connection,'Covariate ');
        condition = str2num(Condition_Split{2}(1));
        
        Connection_Split = strsplit(this_connection,'A{');
        direction = Connection_Split{2}(1);
        if strcmp(direction,'1')
            direction = 'forwards';
        elseif strcmp(direction,'2')
            direction = 'backwards';
        end
        to = str2num(Connection_Split{2}(4));
        from = str2num(Connection_Split{2}(6));
        
        
        if condition == 1
            if BMA_Overall.Ep(these_differences(this_difference),this_contrast)>0
                disp([new_groupstodo{find(template_PEB.M.X(:,this_contrast)==1)} ' stronger than ' new_groupstodo{find(template_PEB.M.X(:,this_contrast)==-1)} ' for ' direction ' ' source_names{from} ' to ' source_names{to} ' with posterior probability ' num2str(BMA_Overall.Pp(these_differences(this_difference),this_contrast))])
            else
                disp([new_groupstodo{find(template_PEB.M.X(:,this_contrast)==-1)} ' stronger than ' new_groupstodo{find(template_PEB.M.X(:,this_contrast)==1)} ' for ' direction ' ' source_names{from} ' to ' source_names{to} ' with posterior probability ' num2str(BMA_Overall.Pp(these_differences(this_difference),this_contrast))])
            end
        else
            disp(['Interaction between ' new_groupstodo{find(template_PEB.M.X(:,this_contrast)==-1)} ' and ' new_groupstodo{find(template_PEB.M.X(:,this_contrast)==1)} ' for ' direction ' ' source_names{from} ' to ' source_names{to} ' with posterior probability ' num2str(BMA_Overall.Pp(these_differences(this_difference),this_contrast))])
        end
        
        %Now repeat combining FTD and AD
        merged_groups = group;
        merged_groups(merged_groups==4) = 2; % combine nfvPPA and bvFTD
        merged_groups(merged_groups==5) = 3; % combine pca and ADMCI
        
        this_fig = figure;
        set(gcf,'Position',[100 100 1600 800]);
        set(gcf, 'PaperPositionMode', 'auto');
        clear linehandle legend_text
        % Reminder: size(all_granger_data{1}) =     8     8     1    20     2   123
        % From, to, time_window, foi, STD-DVT, subj
        for this_subplot = 1:4*length(all_granger_data)
            subplot(4,length(all_granger_data),this_subplot)
            this_measure = mod(this_subplot,length(all_granger_data));
            if this_measure == 0;  this_measure = length(all_granger_data); end
            this_connectivity_contrast = ceil(this_subplot/length(all_granger_data));
            
            if this_connectivity_contrast == 1
                % First look at the overall directionality by group - de-meaned
                % by individual, not meaningful for the icoh contrasts
                difference_demeaned = demeaned_all_granger_data{this_measure}(from,to,:,:,:,:)-demeaned_all_granger_data{this_measure}(to,from,:,:,:,:);
                hold on
                plot([0 ceil(max(foi))],[0 0],'k--','LineWidth',1);
                xlim([0 ceil(max(foi)/10)*10])
                legend_text = cell(1,length(merged_groups));
                for grp = 1:length(unique(merged_groups))
                    if grp == find(template_PEB.M.X(:,this_contrast)==1) || grp == find(template_PEB.M.X(:,this_contrast)==-1)
                        linehandle(grp) = stdshade_TEC_cmap(squeeze(mean(difference_demeaned(1,1,:,:,:,merged_groups==grp),5))',0.2,cmap(grp,:),foi,1,1,':');
                        legend_text{grp} = [new_groupstodo{grp}];
                    end
                end
                if this_measure == length(all_granger_data)
                    legend(linehandle(~cellfun('isempty',legend_text)),legend_text(~cellfun('isempty',legend_text)), 'Interpreter', 'none');
                end
                title(['Directionality ' analysis_type{this_measure}], 'Interpreter', 'none')
                xlabel('Frequency, Hz')
            elseif this_connectivity_contrast == 2
                % Next look at the overall connection strength by group -
                % confounded by SNR for Granger and perhaps pdc
                hold on
                plot([0 ceil(max(foi))],[0 0],'k--','LineWidth',1);
                xlim([0 ceil(max(foi)/10)*10])
                for grp = 1:length(unique(merged_groups))
                    if grp == find(template_PEB.M.X(:,this_contrast)==1) || grp == find(template_PEB.M.X(:,this_contrast)==-1)
                        stdshade_TEC_cmap(squeeze(mean(all_granger_data{this_measure}(from,to,:,:,:,merged_groups==grp),5))',0.2,cmap(grp,:),foi,1,1,':');
                    end
                    if grp>1
                        comp_grp = grp-1;
                    else
                        comp_grp = 3;
                    end
                    for this_band = 1:size(Frequency_bands,1)
                        these_fois = foi>Frequency_bands(this_band,1)&foi<Frequency_bands(this_band,2);
                        [~, pval] = ttest2(squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,merged_groups==comp_grp),5)),1)),squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,merged_groups==grp),5)),1)),'vartype','unequal');
                        if pval <0.05
                            if mean(squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,merged_groups==grp),5)),1)))>mean(squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,merged_groups==comp_grp),5)),1)))
                                disp([new_groupstodo{grp} ' stronger than ' new_groupstodo{comp_grp} ' in the ' Frequence_band_names{this_band} ' band using metric ' analysis_type{this_measure} ' p= ' num2str(pval)])
                            else
                                disp([new_groupstodo{comp_grp} ' stronger than ' new_groupstodo{grp} ' in the ' Frequence_band_names{this_band} ' band using metric ' analysis_type{this_measure} ' p= ' num2str(pval)])
                            end
                        end
                    end
                end
                title(['Strength ' analysis_type{this_measure}], 'Interpreter', 'none')
                xlabel('Frequency, Hz')
            elseif this_connectivity_contrast == 3
                % Next look at the modulation by STD-DVT
                hold on
                plot([0 ceil(max(foi))],[0 0],'k--','LineWidth',1);
                xlim([0 ceil(max(foi)/10)*10])
                for grp = 1:length(unique(merged_groups))
                    if grp == find(template_PEB.M.X(:,this_contrast)==1) || grp == find(template_PEB.M.X(:,this_contrast)==-1)
                        stdshade_TEC_cmap(squeeze(all_mismatch_contrasts{this_measure}(from,to,:,:,merged_groups==grp))',0.2,cmap(grp,:),foi,1,1,':');
                    end
                    if grp>1
                        comp_grp = grp-1;
                    else
                        comp_grp = 3;
                    end
                    for this_band = 1:size(Frequency_bands,1)
                        these_fois = foi>Frequency_bands(this_band,1)&foi<Frequency_bands(this_band,2);
                        [~, pval] = ttest2(squeeze(mean(squeeze(all_mismatch_contrasts{this_measure}(from,to,:,these_fois,merged_groups==comp_grp)),1)),squeeze(mean(squeeze(all_mismatch_contrasts{this_measure}(from,to,:,these_fois,merged_groups==grp)),1)),'vartype','unequal');
                        if pval <0.05
                            disp([new_groupstodo{grp} ' interacts with ' new_groupstodo{comp_grp} ' in the ' Frequence_band_names{this_band} ' band using metric ' analysis_type{this_measure} ' p= ' num2str(pval)])
                        end
                    end
                end
                title(['STD-DVT ' analysis_type{this_measure}], 'Interpreter', 'none')
                xlabel('Frequency, Hz')
            elseif this_connectivity_contrast == 4
                % Finally look at the power for the reverse direction of the same pair -
                % only meaningful for some measures
                % Next look at the overall connection strength by group -
                % confounded by SNR for Granger and perhaps pdc
                hold on
                plot([0 ceil(max(foi))],[0 0],'k--','LineWidth',1);
                xlim([0 ceil(max(foi)/10)*10])
                for grp = 1:length(unique(merged_groups))
                    if grp == find(template_PEB.M.X(:,this_contrast)==1) || grp == find(template_PEB.M.X(:,this_contrast)==-1)
                        stdshade_TEC_cmap(squeeze(mean(all_granger_data{this_measure}(to,from,:,:,:,merged_groups==grp),5))',0.2,cmap(grp,:),foi,1,1,':');
                    end
                end
                title(['Strength ' analysis_type{this_measure}], 'Interpreter', 'none')
                xlabel('Frequency, Hz')
                
            end
        end
        suptitle(['Connectivity from ' source_names{from}  ' to ' source_names{to} ' combined groups']);
        if ceil(max(foi))> 45
            savestring = ['./figures/' source_names{from} '_' source_names{to} '_Multifig_combined_highfreq_focused.pdf'];
        else
            savestring = ['./figures/' source_names{from} '_' source_names{to} '_Multifig_combined_focused.pdf'];
        end
        savestring = strrep(savestring,' ','_');
        if save_figures
            print(savestring,'-depsc','-painters'); %eval(['export_fig ' savestring ' -transparent']);
            eval(['export_fig ' savestring(1:end-3) 'png -transparent']);
        end
        
    end
end

cd([dirname_DCM 'PEB_secondlevel'])
load('PEB_D_Overall_combined.mat')

for this_contrast = 2:size(template_PEB.M.X,2)
    these_differences = find(BMA_Overall.Pp(:,this_contrast)>thresh);
    
    for this_difference = 1:length(these_differences)
        
        delays = {'local','cortico-cortical','cor-thalamo-cor'};
        
        this_connection = BMA_Overall.Pnames{these_differences(this_difference)};
        Condition_Split = strsplit(this_connection,'Covariate ');
        condition = str2num(Condition_Split{2}(1));
        
        Connection_Split = strsplit(this_connection,'D(');
        connection = delays{str2num(Connection_Split{2}(1))};
        
        if condition == 1
            if BMA_Overall.Ep(these_differences(this_difference),this_contrast)>0
                disp([new_groupstodo{find(template_PEB.M.X(:,this_contrast)==1)} ' longer than ' new_groupstodo{find(template_PEB.M.X(:,this_contrast)==-1)} ' for ' connection ' with posterior probability ' num2str(BMA_Overall.Pp(these_differences(this_difference),this_contrast))])
            else
                disp([new_groupstodo{find(template_PEB.M.X(:,this_contrast)==-1)} ' longer than ' new_groupstodo{find(template_PEB.M.X(:,this_contrast)==1)} ' for ' connection ' with posterior probability ' num2str(BMA_Overall.Pp(these_differences(this_difference),this_contrast))])
            end
        else
            disp(['Interaction between ' new_groupstodo{find(template_PEB.M.X(:,this_contrast)==-1)} ' and ' new_groupstodo{find(template_PEB.M.X(:,this_contrast)==1)} ' for ' connection ' with posterior probability ' num2str(BMA_Overall.Pp(these_differences(this_difference),this_contrast))])
        end
        
    end
end

cd(thisdir)

