function tc_batch_SPM(filetype,subjects,pathstem,p)

%% Initialise path and subject definitions

% For info: Trial types
% STD = 1
% DVT = 2
% Loc = 3
% Int = 4
% Dur = 5
% Gap = 6
% Freq = 7
% Loc_L = 8
% Freq_hi = 9
% Int_hi = 10
% Loc_R = 11
% Freq_lo = 12
% Int_lo = 13

%conditions = {'STD','DVT','location','intensity','duration','gap','frequency','location_L','frequency_high','intensity_high','location_R','frequency_low','intensity_low'};
conditions = p.conditions;
ncond = length(conditions);

%% Configure

%filetype = 'PfmraedfffM';
filetypesplit = 'notneededhere';
%modality = {'MEGMAG' 'MEGCOMB'};
modality = p.mod;
imagetype = {'sm_'};
%p.windows = [-100 500];

if ~exist('rejecteeg','var')&&~any(ismember(modality,'EEG'))
    rejecteeg=ones(size(subjects));
elseif exist('rejecteeg','var')&&any(ismember(modality,'EEG'))
    rejecteeg=zeros(size(subjects));
end
rejecteeg = num2cell(rejecteeg);

outputstem = [pathstem '/ERP_stats'];

% Contrasts (don't specify if not needed)
contrasts = {};
cnt = 0;

all_group_combinations = flipud(unique(perms([1, -1, zeros(1,length(p.diagnosis_list)-2)]),'rows','stable'));
all_group_combinations_strings = cell(1,size(all_group_combinations,1));
for this_grp = 1:size(all_group_combinations,1)
    all_group_combinations_strings{this_grp} = sprintf('%s minus %s',p.diagnosis_list{all_group_combinations(this_grp,:)==1},p.diagnosis_list{all_group_combinations(this_grp,:)==-1});
end
% all_condition_combinations = flipud(unique(perms([-1, zeros(1,4)]),'rows','stable'));
% all_condition_combinations = [ones(size(all_condition_combinations,1),1),zeros(size(all_condition_combinations,1),1),all_condition_combinations,zeros(size(all_condition_combinations,1),ncond-2-size(all_condition_combinations,2))];
all_condition_combinations = [ones(length(conditions)-1,1),-eye(length(conditions)-1)];
all_condition_combinations = [all_condition_combinations; -all_condition_combinations];

vs_control_combinations = [ones(length(p.diagnosis_list)-1,1),-eye(length(p.diagnosis_list)-1)];
vs_control_combinations = [vs_control_combinations;-vs_control_combinations];
vs_control_combinations_strings = cell(1,size(vs_control_combinations,1));
for this_grp = 1:size(vs_control_combinations,1)
    vs_control_combinations_strings{this_grp} = sprintf('%s minus %s',p.diagnosis_list{vs_control_combinations(this_grp,:)==1},p.diagnosis_list{vs_control_combinations(this_grp,:)==-1});
end

all_condition_combinations_strings = cell(1,size(all_condition_combinations,1));
for this_cond = 1:size(all_condition_combinations,1)
    all_condition_combinations_strings{this_cond} = sprintf('%s minus %s',conditions{all_condition_combinations(this_cond,:)==1},conditions{all_condition_combinations(this_cond,:)==-1});
end


%% Contrasts (Combined SPM for patients/controls)

if length(conditions)==2
cnt = cnt + 1;
contrasts{cnt}.name = 'Mismatch response (All)';
contrasts{cnt}.c = kron(ones(1,length(p.diagnosis_list)),[1,-1,zeros(1,ncond-2)]);
contrasts{cnt}.type = 'F';

cnt = cnt + 1;
contrasts{cnt}.name = ['Mismatch response controls minus allpats']; %NB weights each patient group equally regardless of size
contrasts{cnt}.c = kron([1, -(ones(1,length(p.diagnosis_list)-1)/(length(p.diagnosis_list)-1))],[1,-1,zeros(1,ncond-2)]);
contrasts{cnt}.type = 'T';

cnt = cnt + 1;
contrasts{cnt}.name = ['Mismatch response allpats minus controls'];
contrasts{cnt}.c = kron([-1, (ones(1,length(p.diagnosis_list)-1)/(length(p.diagnosis_list)-1))],[1,-1,zeros(1,ncond-2)]);
contrasts{cnt}.type = 'T';

cnt = cnt + 1;
contrasts{cnt}.name = ['Mismatch response controls vs eachgroup'];
contrasts{cnt}.c = kron(vs_control_combinations(1:size(vs_control_combinations,1)/2,:),[1,-1,zeros(1,ncond-2)]);
contrasts{cnt}.type = 'F';
end

for this_grp = 1:size(all_group_combinations,1)
    for this_cond = 1:size(all_condition_combinations,1)
        cnt = cnt + 1;
        contrasts{cnt}.name = [all_condition_combinations_strings{this_cond} '; ' all_group_combinations_strings{this_grp}];
        contrasts{cnt}.c = kron(all_group_combinations(this_grp,:),all_condition_combinations(this_cond,:));
        contrasts{cnt}.type = 'T';
    end
end

% 
% for this_grp = 1:size(vs_control_combinations,1)
%     cnt = cnt + 1;
%     contrasts{cnt}.name = ['Mismatch response ' vs_control_combinations_strings{this_grp}];
%     contrasts{cnt}.c = kron(vs_control_combinations(this_grp,:),[1,-1,zeros(1,ncond-2)]);
%     contrasts{cnt}.type = 'T';
% end
% 
% for this_cond = 1:size(all_condition_combinations,1)
%     cnt = cnt + 1;
%     contrasts{cnt}.name = [all_condition_combinations_strings{this_cond} ' controls vs eachgroup'];
%     contrasts{cnt}.c = kron(vs_control_combinations(1:size(vs_control_combinations,1)/2,:),all_condition_combinations(this_cond,:));
%     contrasts{cnt}.type = 'F';
% end
% 
% for this_grp = 1:size(vs_control_combinations,1)
%     for this_cond = 1:size(all_condition_combinations,1)
%         cnt = cnt + 1;
%         contrasts{cnt}.name = [all_condition_combinations_strings{this_cond} '; ' vs_control_combinations_strings{this_grp}];
%         contrasts{cnt}.c = kron(vs_control_combinations(this_grp,:),all_condition_combinations(this_cond,:));
%         contrasts{cnt}.type = 'T';
%     end
% end


%% Estimate models


%for img=1:length(imagetype)
img = 1;
%for wind = 1:length(p.windows)

for wind = 1
    for m=1:length(modality)
        files = {};
        %for m = 3
        for s=1:length(subjects)
            this_subj_id = strsplit(subjects{s},'/');
            this_subj_id = this_subj_id{2};

            for c=1:length(conditions)
                if strcmp(modality{m},'EEG')
                    if rejecteeg{s} == 1
                        %files{p.group(s)}{s}{c} = [];
                    else
                        files{p.group(s)}{s}{c} = strjoin([pathstem subjects{s} '/' modality{m} filetype this_subj_id '/' imagetype 'condition_' conditions{c} '.nii'],'');
                        if exist(files{p.group(s)}{s}{c},'file')
                        else
                            error([files{p.group(s)}{s}{c} ' does not exist'])
                            %files{p.group(s)}{s}{c} = strjoin([pathstem subjects{s} '/' modality{m} filetypesplit '/' imagetype 'condition_' conditions{c} '.nii'],'');
                        end
                    end
                    
                else
                    files{p.group(s)}{s}{c} = strjoin([pathstem subjects{s} '/' modality{m} filetype this_subj_id '/' imagetype 'condition_' conditions{c} '.nii'],'');
                    if exist(files{p.group(s)}{s}{c},'file')
                    else
                        error([files{p.group(s)}{s}{c} ' does not exist'])
                        %files{p.group(s)}{s}{c} = strjoin([pathstem subjects{s} '/' modality{m} filetypesplit '/' imagetype 'condition_' conditions{c} '.nii'],'');
                    end
                end
            end
            
            
        end
        
        for this_grp = 1:length(p.diagnosis_list)
            files{this_grp} = files{this_grp}(~cellfun(@isempty,files{this_grp}));
        end
        
        % set up input structure for batch_spm_anova_vES
        S.imgfiles = files;
        if ~isempty(strfind(S.imgfiles{1}{1}{1},'tf_'))
            outputfullpath = [outputstem imagetype{img} '/TF_combined_' strjoin(p.conditions,'_') num2str(p.windows(wind,1)) '_' num2str(p.windows(wind,2)) '_' modality{m}];
            mskname = [pathstem modality{m} sprintf(['_TF_mask_%d_%dms.img'],p.windows(wind,1),p.windows(wind,2))];
        else
            outputfullpath = [outputstem imagetype{img} '/combined_' strjoin(p.conditions,'_') num2str(p.windows(wind,1)) '_' num2str(p.windows(wind,2)) '_' modality{m}];
            mskname = [pathstem modality{m} sprintf(['_mask_%d_%dms.img'],p.windows(wind,1),p.windows(wind,2))];
        end
        if ~exist(outputfullpath)
            mkdir(outputfullpath);
        end
        S.outdir = outputfullpath;
        S.uUFp = 1; % for M/EEG only
        %S.nsph_flag = 0;
%         if strncmp(modality{m},'time_',5)
%             %mskname = [pathstem modality{m}(6:end)
%             %'_1D_mask_0_800ms.img']; No need for mask - images created
%             %with restricted time window
%         else
%             
%             %mskname = [pathstem modality{m} '_mask_-100_800ms.img'];
%         end
        if exist('mskname'); S.maskimg = mskname; end;
        if exist('contrasts'); S.contrasts = contrasts; end;
        if exist('covariates'); S.user_regs = covariates; end;
        
        % estimate model and compute contrasts
        batch_spm_anova_es(S);
        %pause
    end
end

%end
