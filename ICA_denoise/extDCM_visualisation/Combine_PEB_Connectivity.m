function Combine_PEB_Connectivity(dirname_DCM,diagnosis_list,source_names,thresh,Participant)
%A script for plotting the results of extDCM across all diagnoses by
%inter-regional connection
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
timewins = [0 500];
%topfreqband = 49;
topfreqband = 20;

closeafter = 1;

analysis_type = {};
analysis_type{end+1} = 'Granger';
analysis_type{end+1} = 'icoh';
analysis_type{end+1} = 'plv';
analysis_type{end+1} = 'pdc';
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
            case {'icoh', 'partial_icoh'}
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
                disp([diagnosis_list{find(template_PEB.M.X(:,this_contrast)==1)} ' stronger than ' diagnosis_list{find(template_PEB.M.X(:,this_contrast)==-1)} ' for ' direction ' ' source_names{from} ' to ' source_names{to}])
            else
                disp([diagnosis_list{find(template_PEB.M.X(:,this_contrast)==-1)} ' stronger than ' diagnosis_list{find(template_PEB.M.X(:,this_contrast)==1)} ' for ' direction ' ' source_names{from} ' to ' source_names{to}])
            end
        else
            disp(['Interaction between ' diagnosis_list{find(template_PEB.M.X(:,this_contrast)==-1)} ' and ' diagnosis_list{find(template_PEB.M.X(:,this_contrast)==1)} ' for ' direction ' ' source_names{from} ' to ' source_names{to}])
        end
        
        this_fig = figure;
        set(gcf,'Position',[100 100 1600 800]);
        set(gcf, 'PaperPositionMode', 'auto');
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
                plot([0 40],[0 0],'k--','LineWidth',1);
                hold on
                legend_text = cell(1,2*length(groupstodo));
                for grp = 1:2:2*length(groupstodo)
                    linehandle(grp) = stdshade_TEC_cmap(squeeze(mean(difference_demeaned(1,1,:,:,:,group==ceil(grp/2)),5))',0.2,cmap(ceil(grp/2),:),foi,1,1,':');
                    legend_text{grp} = [groupstodo{ceil(grp/2)}];
                end
                if this_measure == length(all_granger_data)
                    legend(linehandle(~cellfun('isempty',legend_text)),legend_text(~cellfun('isempty',legend_text)));
                end
                title(['Directionality ' analysis_type{this_measure}], 'Interpreter', 'none')
                xlabel('Frequency, Hz')
            elseif this_connectivity_contrast == 2
                % Next look at the overall connection strength by group -
                % confounded by SNR for Granger and perhaps pdc
            end
        end
                
        
    end
end

load('PEB_D_Overall.mat')

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
                disp([diagnosis_list{find(template_PEB.M.X(:,this_contrast)==1)} ' longer than ' diagnosis_list{find(template_PEB.M.X(:,this_contrast)==-1)} ' for ' connection])
            else
                disp([diagnosis_list{find(template_PEB.M.X(:,this_contrast)==-1)} ' longer than ' diagnosis_list{find(template_PEB.M.X(:,this_contrast)==1)} ' for ' connection])
            end
        else
            disp(['Interaction between ' diagnosis_list{find(template_PEB.M.X(:,this_contrast)==-1)} ' and ' diagnosis_list{find(template_PEB.M.X(:,this_contrast)==1)} ' for ' connection])
        end
        
    end
end

cd(thisdir)
