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
    
end

cd(thisdir)

group = all_subjs(:,1)';


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
