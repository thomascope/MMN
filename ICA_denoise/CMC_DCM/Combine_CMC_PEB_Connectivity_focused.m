function Combine_CMC_PEB_Connectivity_focused(dirname_DCM,diagnosis_list,source_names,thresh,Participant,PEB_focus)
%A script for plotting the results of extDCM across all diagnoses by
%inter-regional connection
save_figures = 0;
do_nonsigs = 0;
defaultStream=RandStream('mt19937ar','Seed',15); % For later permutation tests - ensure that the same seed is used for reproducability

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
cd([dirname_DCM 'PEB_secondlevel'])

load(['PEB_' PEB_focus '_' cat(2,diagnosis_list{:}) '.mat'])
template_PEB = PEB;

assert(all(template_PEB.M.X(:,1)==1),'The first column of the PEB of PEBs contrast should be all ones, check please.')

base_datapathstem = '/imaging/tc02/Holly_MMN/Coherence_Connectivity_Integrated_LOR/crosshem/'; %With sLORETA source reconstruction
addpath(['/group/language/data/thomascope/vespa/SPM12version/Standalone preprocessing pipeline/tc_source_stats/ojwoodford-export_fig-216b30e']);
addpath(['/group/language/data/thomascope/MMN/ICA_denoise/stdshade']);
addpath('/imaging/tc02/toolboxes/rsatoolbox/Engines/');

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
start_times = 0;
end_times = 500;
topfreqband = 40; % Include gamma
%topfreqband = 20;

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
main_effect_significances = [];
interaction_significances = [];
permutation_main_effect_significances = [];
permutation_interaction_significances = [];

for this_contrast = 2:size(template_PEB.M.X,2)
    these_differences = find(BMA.Pp(:,this_contrast)>thresh);
    
    for this_difference = 1:length(these_differences)
        
        this_connection = BMA.Pnames{these_differences(this_difference)};
        %         Condition_Split = strsplit(this_connection,'Covariate ');
        %         condition = str2num(Condition_Split{2}(1));
        
        if strcmp(PEB_focus,'A')
            Connection_Split = strsplit(this_connection,[PEB_focus '{']);
            if str2num(Connection_Split{2}(1)) > 2
                condition = str2num(Connection_Split{2}(1));
                direction = '2';
            else
                condition = str2num(Connection_Split{2}(1));
                direction = '1';
            end
            if strcmp(direction,'1')
                direction = 'forwards';
            elseif strcmp(direction,'2')
                direction = 'backwards';
            end
            to = str2num(Connection_Split{2}(4));
            from = str2num(Connection_Split{2}(6));
            
            if BMA.Ep(these_differences(this_difference),this_contrast)>0
                disp([Participant{min(find(template_PEB.M.X(:,this_contrast)==1))}.diag ' stronger than ' Participant{min(find(template_PEB.M.X(:,this_contrast)==-1))}.diag ' for ' direction ' ' source_names{from} ' to ' source_names{to} ' ' cell_pops{condition} ' with posterior probability ' num2str(BMA.Pp(these_differences(this_difference),this_contrast))])
            else
                disp([Participant{min(find(template_PEB.M.X(:,this_contrast)==-1))}.diag ' stronger than ' Participant{min(find(template_PEB.M.X(:,this_contrast)==1))}.diag ' for ' direction ' ' source_names{from} ' to ' source_names{to} ' ' cell_pops{condition} ' with posterior probability ' num2str(BMA.Pp(these_differences(this_difference),this_contrast))])
            end
            
        elseif strcmp(PEB_focus,'B')
            Connection_Split = strsplit(this_connection,[PEB_focus '{']);
            to = str2num(Connection_Split{2}(4));
            from = str2num(Connection_Split{2}(6));
            
            if BMA.Ep(these_differences(this_difference),this_contrast)>0
                disp([Participant{min(find(template_PEB.M.X(:,this_contrast)==1))}.diag ' stronger than ' Participant{min(find(template_PEB.M.X(:,this_contrast)==-1))}.diag ' for modulation of ' source_names{from} ' to ' source_names{to} ' with posterior probability ' num2str(BMA.Pp(these_differences(this_difference),this_contrast))])
            else
                disp([Participant{min(find(template_PEB.M.X(:,this_contrast)==-1))}.diag ' stronger than ' Participant{min(find(template_PEB.M.X(:,this_contrast)==1))}.diag ' for modulation of ' source_names{from} ' to ' source_names{to} ' with posterior probability ' num2str(BMA.Pp(these_differences(this_difference),this_contrast))])
            end
            
        else
            error('So far only implemented for spm_fx_CMC A and B matrices')
        end
        
        if from == to
            disp('Intrinsic connection, moving on')
            continue
        end
        
        this_fig = figure;
        set(gcf,'Position',[100 100 1600 800]);
        set(gcf, 'PaperPositionMode', 'auto');
        clear linehandle legend_text
        % Reminder: size(all_granger_data{1}) =     8     8     1    20     2   123
        % From, to, time_window, foi, STD-DVT, subj
        if strcmp(PEB_focus,'A')
            main_effect_significances(end+1) = 0;
            permutation_main_effect_significances(end+1) = 0;
        else
            interaction_significances(end+1) = 0;
            permutation_interaction_significances(end+1) = 0;
        end
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
                    if strcmp(diagnosis_list{ceil(grp/2)},Participant{min(find(template_PEB.M.X(:,this_contrast)==1))}.diag) || strcmp(diagnosis_list{ceil(grp/2)},Participant{min(find(template_PEB.M.X(:,this_contrast)==-1))}.diag)
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
                    if strcmp(diagnosis_list{ceil(grp/2)},Participant{min(find(template_PEB.M.X(:,this_contrast)==1))}.diag) || strcmp(diagnosis_list{ceil(grp/2)},Participant{min(find(template_PEB.M.X(:,this_contrast)==-1))}.diag)
                        stdshade_TEC_cmap(squeeze(mean(all_granger_data{this_measure}(from,to,:,:,:,group==ceil(grp/2)),5))',0.2,cmap(ceil(grp/2),:),foi,1,1,':');
                    end
                end
                if strcmp(PEB_focus,'A')
                    data_for_test = [];
                    data_for_test_2 = [];
                    %First standard parametric tests per frequency band - works well but perhaps too statistically lenient
                    %Then permutation tests - more robust!
                    for this_band = 1:size(Frequency_bands,1)
                        these_fois = foi>Frequency_bands(this_band,1)&foi<Frequency_bands(this_band,2);
                        [~, pval] = ttest2(squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==1)),5)),1)),squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==-1)),5)),1)),'vartype','unequal');
                        if isnan(pval)
                            pause
                        end
                        permutation_p = permutationTest(squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==1)),5)),1)),squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==-1)),5)),1)),num_perms);
                        if pval < post_thresh || permutation_p < post_thresh
                            if pval < post_thresh
                                main_effect_significances(end) = 1;
                            end
                            if permutation_p < post_thresh
                                permutation_main_effect_significances(end) = 1;
                            end
                            if mean(squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==1)),5)),1)))>mean(squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==-1)),5)),1)))
                                disp([Participant{min(find(template_PEB.M.X(:,this_contrast)==1))}.diag ' stronger than ' Participant{min(find(template_PEB.M.X(:,this_contrast)==-1))}.diag ' in the ' Frequence_band_names{this_band} ' band using metric ' analysis_type{this_measure} ' parametric p= ' num2str(pval) ' permutation p= ' num2str(permutation_p)])
                            else
                                disp([Participant{min(find(template_PEB.M.X(:,this_contrast)==-1))}.diag ' stronger than ' Participant{min(find(template_PEB.M.X(:,this_contrast)==1))}.diag ' in the ' Frequence_band_names{this_band} ' band using metric ' analysis_type{this_measure} ' parametric p= ' num2str(pval) ' permutation p= ' num2str(permutation_p)])
                            end
                        end
                        data_for_test = [data_for_test;squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==1)),5)),1))];
                        data_for_test_2 = [data_for_test_2;squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==-1)),5)),1))];
                    end
                    %                     %Next at an attempt at permutation tests across frequency bands - I don't think that this works
                    %                     %properly so I have left the code here but changed the number of permutations to
                    %                     %10 in order to make it run quickly. Change back to 5000 or so to implement this method if you trust it.
                    %                     num_perms_t2 = 10;
                    %                     pval=mult_comp_perm_t2(data_for_test',data_for_test_2',num_perms_t2,0,0.05,0,'w',0,defaultStream.State);
                    %                     for this_band = 1:size(Frequency_bands,1)
                    %                         if pval(this_band) < post_thresh
                    %                             permutation_main_effect_significances(end) = 1;
                    %                             if mean(squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==1)),5)),1)))>mean(squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==-1)),5)),1)))
                    %                                 disp([Participant{min(find(template_PEB.M.X(:,this_contrast)==1))}.diag ' stronger than ' Participant{min(find(template_PEB.M.X(:,this_contrast)==-1))}.diag ' in the ' Frequence_band_names{this_band} ' band using metric ' analysis_type{this_measure} ' permutation p= ' num2str(pval(this_band))])
                    %                             else
                    %                                 disp([Participant{min(find(template_PEB.M.X(:,this_contrast)==-1))}.diag ' stronger than ' Participant{min(find(template_PEB.M.X(:,this_contrast)==1))}.diag ' in the ' Frequence_band_names{this_band} ' band using metric ' analysis_type{this_measure} ' permutation p= ' num2str(pval(this_band))])
                    %                             end
                    %                         end
                    %                     end
                    
                end
                title(['Strength ' analysis_type{this_measure}], 'Interpreter', 'none')
                xlabel('Frequency, Hz')
            elseif this_connectivity_contrast == 3
                % Next look at the modulation by STD-DVT
                hold on
                plot([0 ceil(max(foi))],[0 0],'k--','LineWidth',1);
                xlim([0 ceil(max(foi)/10)*10])
                for grp = 1:2:2*length(groupstodo)
                    if strcmp(diagnosis_list{ceil(grp/2)},Participant{min(find(template_PEB.M.X(:,this_contrast)==1))}.diag) || strcmp(diagnosis_list{ceil(grp/2)},Participant{min(find(template_PEB.M.X(:,this_contrast)==-1))}.diag)
                        stdshade_TEC_cmap(squeeze(all_mismatch_contrasts{this_measure}(from,to,:,:,group==ceil(grp/2)))',0.2,cmap(ceil(grp/2),:),foi,1,1,':');
                    end
                end
                if strcmp(PEB_focus,'B')
                    data_for_test = [];
                    data_for_test_2 = [];
                    for this_band = 1:size(Frequency_bands,1)
                        these_fois = foi>Frequency_bands(this_band,1)&foi<Frequency_bands(this_band,2);
                        [~, pval] = ttest2(squeeze(mean(squeeze(all_mismatch_contrasts{this_measure}(from,to,:,these_fois,find(template_PEB.M.X(:,this_contrast)==1))),1)),squeeze(mean(squeeze(all_mismatch_contrasts{this_measure}(from,to,:,these_fois,find(template_PEB.M.X(:,this_contrast)==-1))),1)),'vartype','unequal');
                        permutation_p = permutationTest(squeeze(mean(squeeze(all_mismatch_contrasts{this_measure}(from,to,:,these_fois,find(template_PEB.M.X(:,this_contrast)==1))),1)),squeeze(mean(squeeze(all_mismatch_contrasts{this_measure}(from,to,:,these_fois,find(template_PEB.M.X(:,this_contrast)==-1))),1)),num_perms);
                        if pval < post_thresh || permutation_p < post_thresh
                            if pval < post_thresh
                                interaction_significances(end) = 1;
                            end
                            if permutation_p < post_thresh
                                permutation_interaction_significances(end) = 1;
                            end
                            disp([Participant{min(find(template_PEB.M.X(:,this_contrast)==1))}.diag ' interacts with ' Participant{min(find(template_PEB.M.X(:,this_contrast)==-1))}.diag ' in the ' Frequence_band_names{this_band} ' band using metric ' analysis_type{this_measure} ' parametric p= ' num2str(pval) ' permutation p= ' num2str(permutation_p)])
                        end
                        data_for_test = [data_for_test;squeeze(mean(squeeze(all_mismatch_contrasts{this_measure}(from,to,:,these_fois,find(template_PEB.M.X(:,this_contrast)==1))),1))];
                        data_for_test_2 = [data_for_test_2;squeeze(mean(squeeze(all_mismatch_contrasts{this_measure}(from,to,:,these_fois,find(template_PEB.M.X(:,this_contrast)==-1))),1))];
                    end
                    
                    %                     pval=mult_comp_perm_t2(data_for_test',data_for_test_2',num_perms_t2,0,0.05,0,'w',0,defaultStream.State);
                    %                     for this_band = 1:size(Frequency_bands,1)
                    %                         if pval(this_band) < post_thresh
                    %                             permutation_interaction_significances(end) = 1;
                    %                             disp([Participant{min(find(template_PEB.M.X(:,this_contrast)==1))}.diag ' interacts with ' Participant{min(find(template_PEB.M.X(:,this_contrast)==-1))}.diag ' in the ' Frequence_band_names{this_band} ' band using metric ' analysis_type{this_measure} ' permutation p= ' num2str(pval(this_band))])
                    %                         end
                    %                     end
                    
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
                    if strcmp(diagnosis_list{ceil(grp/2)},Participant{min(find(template_PEB.M.X(:,this_contrast)==1))}.diag) || strcmp(diagnosis_list{ceil(grp/2)},Participant{min(find(template_PEB.M.X(:,this_contrast)==-1))}.diag)
                        stdshade_TEC_cmap(squeeze(mean(all_granger_data{this_measure}(to,from,:,:,:,group==ceil(grp/2)),5))',0.2,cmap(ceil(grp/2),:),foi,1,1,':');
                    end
                end
                title(['Strength ' analysis_type{this_measure}], 'Interpreter', 'none')
                xlabel('Frequency, Hz')
                
            end
        end
        suptitle(['Connectivity from ' source_names{from}  ' to ' source_names{to}]);
        if ceil(max(foi)) > 45
            savestring = ['./CMC_figures/' all_folders{end} '_' cat(2,diagnosis_list{:}) '_' source_names{from} '_' source_names{to} '_Multifig_highfreq_focused.pdf'];
        else
            savestring = ['./CMC_figures/' all_folders{end} '_' cat(2,diagnosis_list{:}) '_' source_names{from} '_' source_names{to} '_Multifig_focused.pdf'];
        end
        
        savestring = strrep(savestring,' ','_');
        if save_figures
            print(savestring,'-depsc','-painters'); %eval(['export_fig ' savestring ' -transparent']);
            eval(['export_fig ' savestring(1:end-3) 'png -transparent']);
        end
        
    end
end

if do_nonsigs %Not yet implemented for CMC - run the Combine_PEB_Connectivity_focused script instead
    overall_main_effect_significances = [];
    overall_interaction_significances = [];
    overall_permutation_main_effect_significances = [];
    overall_permutation_interaction_significances = [];
    % Now find out what proportion of non-DCM PEB connections are significant in coherence/Granger/plv
    for this_contrast = 2:size(template_PEB.M.X,2)
        
        for this_difference = 1:length(BMA.Pnames)
            this_connection = BMA.Pnames{this_difference};
            Condition_Split = strsplit(this_connection,'Covariate ');
            condition = str2num(Condition_Split{2}(1));
            
            Connection_Split = strsplit(this_connection,[PEB_focus '{']);
            direction = Connection_Split{2}(1);
            if strcmp(direction,'1')
                direction = 'forwards';
            elseif strcmp(direction,'2')
                direction = 'backwards';
            end
            to = str2num(Connection_Split{2}(4));
            from = str2num(Connection_Split{2}(6));
            
            if strcmp(PEB_focus,'A')
                overall_main_effect_significances(end+1) = 0;
                overall_permutation_main_effect_significances(end+1) = 0;
                for this_measure = 1:length(all_granger_data)
                    data_for_test = [];
                    data_for_test_2 = [];
                    for this_band = 1:size(Frequency_bands,1)
                        these_fois = foi>Frequency_bands(this_band,1)&foi<Frequency_bands(this_band,2);
                        [~, pval] = ttest2(squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==1)),5)),1)),squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==-1)),5)),1)),'vartype','unequal');
                        permutation_p = permutationTest(squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==1)),5)),1)),squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==-1)),5)),1)),num_perms);
                        if pval < post_thresh || permutation_p < post_thresh
                            if pval < post_thresh
                                overall_main_effect_significances(end) = 1;
                            end
                            if permutation_p < post_thresh
                                overall_permutation_main_effect_significances(end) = 1;
                            end
                            if mean(squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==1)),5)),1)))>mean(squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==-1)),5)),1)))
                                %disp([Participant{min(find(template_PEB.M.X(:,this_contrast)==1))}.diag ' stronger than ' Participant{min(find(template_PEB.M.X(:,this_contrast)==-1))}.diag ' in the ' Frequence_band_names{this_band} ' band using metric ' analysis_type{this_measure} ' p= ' num2str(pval)])
                            else
                                %disp([Participant{min(find(template_PEB.M.X(:,this_contrast)==-1))}.diag ' stronger than ' Participant{min(find(template_PEB.M.X(:,this_contrast)==1))}.diag ' in the ' Frequence_band_names{this_band} ' band using metric ' analysis_type{this_measure} ' p= ' num2str(pval)])
                            end
                        end
                        data_for_test = [data_for_test;squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==1)),5)),1))];
                        data_for_test_2 = [data_for_test_2;squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==-1)),5)),1))];
                    end
                    
                    %                 pval=mult_comp_perm_t2(data_for_test',data_for_test_2',num_perms_t2,0,0.05,0,'w',0,defaultStream.State);
                    %                 for this_band = 1:size(Frequency_bands,1)
                    %                     if pval(this_band) < post_thresh
                    %                         overall_permutation_main_effect_significances(end) = 1;
                    %                         if mean(squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==1)),5)),1)))>mean(squeeze(mean(squeeze(mean(all_granger_data{this_measure}(from,to,:,these_fois,:,find(template_PEB.M.X(:,this_contrast)==-1)),5)),1)))
                    %                             disp([Participant{min(find(template_PEB.M.X(:,this_contrast)==1))}.diag ' stronger than ' Participant{min(find(template_PEB.M.X(:,this_contrast)==-1))}.diag ' in the ' Frequence_band_names{this_band} ' band using metric ' analysis_type{this_measure} ' permutation p= ' num2str(pval(this_band))])
                    %                         else
                    %                             disp([Participant{min(find(template_PEB.M.X(:,this_contrast)==-1))}.diag ' stronger than ' Participant{min(find(template_PEB.M.X(:,this_contrast)==1))}.diag ' in the ' Frequence_band_names{this_band} ' band using metric ' analysis_type{this_measure} ' permutation p= ' num2str(pval(this_band))])
                    %                         end
                    %                     end
                    %
                    %                 end
                end
            elseif strcmp(PEB_focus,'B')
                overall_interaction_significances(end+1) = 0;
                overall_permutation_interaction_significances(end+1) = 0;
                for this_measure = 1:length(all_granger_data)
                    data_for_test = [];
                    data_for_test_2 = [];
                    for this_band = 1:size(Frequency_bands,1)
                        these_fois = foi>Frequency_bands(this_band,1)&foi<Frequency_bands(this_band,2);
                        [~, pval] = ttest2(squeeze(mean(squeeze(all_mismatch_contrasts{this_measure}(from,to,:,these_fois,find(template_PEB.M.X(:,this_contrast)==1))),1)),squeeze(mean(squeeze(all_mismatch_contrasts{this_measure}(from,to,:,these_fois,find(template_PEB.M.X(:,this_contrast)==-1))),1)),'vartype','unequal');
                        permutation_p = permutationTest(squeeze(mean(squeeze(all_mismatch_contrasts{this_measure}(from,to,:,these_fois,find(template_PEB.M.X(:,this_contrast)==1))),1)),squeeze(mean(squeeze(all_mismatch_contrasts{this_measure}(from,to,:,these_fois,find(template_PEB.M.X(:,this_contrast)==-1))),1)),num_perms);
                        if pval < post_thresh || permutation_p < post_thresh
                            if pval < post_thresh
                                overall_interaction_significances(end) = 1;
                            end
                            if permutation_p < post_thresh
                                overall_permutation_interaction_significances(end) = 1;
                            end
                            %disp([Participant{min(find(template_PEB.M.X(:,this_contrast)==1))}.diag ' interacts with ' Participant{min(find(template_PEB.M.X(:,this_contrast)==-1))}.diag ' in the ' Frequence_band_names{this_band} ' band using metric ' analysis_type{this_measure} ' p= ' num2str(pval)])
                        end
                        data_for_test = [data_for_test;squeeze(mean(squeeze(all_mismatch_contrasts{this_measure}(from,to,:,these_fois,find(template_PEB.M.X(:,this_contrast)==1))),1))];
                        data_for_test_2 = [data_for_test_2;squeeze(mean(squeeze(all_mismatch_contrasts{this_measure}(from,to,:,these_fois,find(template_PEB.M.X(:,this_contrast)==-1))),1))];
                    end
                    
                    %                 pval=mult_comp_perm_t2(data_for_test',data_for_test_2',num_perms_t2,0,0.05,0,'w',0,defaultStream.State);
                    %                 for this_band = 1:size(Frequency_bands,1)
                    %                     if pval(this_band) < post_thresh
                    %                         overall_permutation_interaction_significances(end) = 1;
                    %                         disp([Participant{min(find(template_PEB.M.X(:,this_contrast)==1))}.diag ' interacts with ' Participant{min(find(template_PEB.M.X(:,this_contrast)==-1))}.diag ' in the ' Frequence_band_names{this_band} ' band using metric ' analysis_type{this_measure} ' permutation p= ' num2str(pval(this_band))])
                    %                     end
                    %                 end
                    %
                end
            end
        end
    end
end

prop_verified = sum(main_effect_significances) / length(main_effect_significances)
prop_inter_verified = sum(interaction_significances) / length(interaction_significances)
prop_permutation_verified = sum(permutation_main_effect_significances) / length(permutation_main_effect_significances)
prop_permutation_inter_verified = sum(permutation_interaction_significances) / length(permutation_interaction_significances)
% prop_grpcombined_verified = sum(grpcombined_main_effect_significances) / length(grpcombined_main_effect_significances)
% prop_grpcombined_inter_verified = sum(grpcombined_interaction_significances) / length(grpcombined_interaction_significances)
% prop_grpcombined_permutation_verified = sum(grpcombined_permutation_main_effect_significances) / length(grpcombined_permutation_main_effect_significances)
% prop_grpcombined_permutation_inter_verified = sum(grpcombined_permutation_interaction_significances) / length(grpcombined_permutation_interaction_significances)

if do_nonsigs
    prop_overall = (sum(overall_main_effect_significances)-sum(main_effect_significances))/(length(overall_main_effect_significances)-length(main_effect_significances))
    prop_inter_overall = (sum(overall_interaction_significances)-sum(interaction_significances))/(length(overall_interaction_significances)-length(interaction_significances))
    prop_permutation_overall = (sum(overall_permutation_main_effect_significances)-sum(permutation_main_effect_significances))/(length(overall_permutation_main_effect_significances)-length(permutation_main_effect_significances))
    prop_permutation_inter_overall = (sum(overall_permutation_interaction_significances)-sum(permutation_interaction_significances))/(length(overall_permutation_interaction_significances)-length(permutation_interaction_significances))
end

cd(thisdir)

workspaceVars = who;
findVars = strncmp(workspaceVars, 'prop_',5);
indexVars = find(findVars);

save(['./CMC_figures/proportions_' cat(2,diagnosis_list{:}) '_' all_folders{end} '.mat'],workspaceVars{indexVars});