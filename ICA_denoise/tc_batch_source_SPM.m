function tc_batch_source_SPM(filetype,Participant,pathstem,p,inv_meth)
% A function for doing a source space SPM for the MMN data


%% Initialise path and subject definitions
addpath('/group/language/data/thomascope/vespa/SPM12version/Standalone preprocessing pipeline/tc_source_stats');
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
mskname = '/imaging/local/spm/spm8/apriori/grey.nii'; % set to [] if not needed

%% Check the input files exist
image_file_pattern = cell(length(p.time_wind_path),length(Participant),length(p.these_conditions));
image_files = cell(length(p.time_wind_path),length(Participant),length(p.these_conditions));
for this_wind = 1:length(p.time_wind_path)
    for this_subj = 1:length(Participant)
        for this_condition = 1:length(p.these_conditions)
            
            image_file_pattern{this_wind,this_subj,this_condition} = [pathstem filesep Participant{this_subj}.groupfolder filesep Participant{this_subj}.namepostmerge filesep 's_' p.time_wind_path{this_wind} '_' inv_meth '_' filetype Participant{this_subj}.namepostmerge '*_' num2str(p.these_conditions(this_condition)) '.nii'] ;
            
            try
                image_files{this_wind,this_subj,this_condition} = ls(image_file_pattern{this_wind,this_subj,this_condition});
            catch
                error(['File does not exist for ' image_file_pattern{this_wind,this_subj,this_condition}]);
            end
        end
    end
end


%conditions = {'STD','DVT','location','intensity','duration','gap','frequency','location_L','frequency_high','intensity_high','location_R','frequency_low','intensity_low'};
conditions = p.conditions(p.these_conditions);
ncond = length(conditions);

%% Configure

outputstem = [pathstem '/source_stats'];


all_group_combinations = flipud(unique(perms([1, -1, zeros(1,length(p.diagnosis_list)-2)]),'rows','stable'));
all_group_combinations_strings = cell(1,size(all_group_combinations,1));
for this_grp = 1:size(all_group_combinations,1)
    all_group_combinations_strings{this_grp} = sprintf('%s minus %s',p.diagnosis_list{all_group_combinations(this_grp,:)==1},p.diagnosis_list{all_group_combinations(this_grp,:)==-1});
end
% all_condition_combinations = flipud(unique(perms([-1, zeros(1,4)]),'rows','stable'));
% all_condition_combinations = [ones(size(all_condition_combinations,1),1),zeros(size(all_condition_combinations,1),1),all_condition_combinations,zeros(size(all_condition_combinations,1),ncond-2-size(all_condition_combinations,2))];
all_condition_combinations = [ones(length(conditions)-1,1),-eye(length(conditions)-1)];
all_condition_combinations = [all_condition_combinations; -all_condition_combinations];

all_F_condition_contrasts = unique(perms([1,-1, zeros(1,length(conditions)-2)]),'rows');
all_F_condition_contrasts = -all_F_condition_contrasts(1:end/2,:);

vs_control_combinations = [ones(length(p.diagnosis_list)-1,1),-eye(length(p.diagnosis_list)-1)];
vs_control_combinations = [vs_control_combinations;-vs_control_combinations];
vs_control_combinations_strings = cell(1,size(vs_control_combinations,1));
for this_grp = 1:size(vs_control_combinations,1)
    vs_control_combinations_strings{this_grp} = sprintf('%s minus %s',p.diagnosis_list{vs_control_combinations(this_grp,:)==1},p.diagnosis_list{vs_control_combinations(this_grp,:)==-1});
end
vs_control_combinations_strings{end+1} = 'All participants';
vs_control_combinations = [vs_control_combinations;ones(1,length(p.diagnosis_list))];

all_condition_combinations_strings = cell(1,size(all_condition_combinations,1));
for this_cond = 1:size(all_condition_combinations,1)
    all_condition_combinations_strings{this_cond} = sprintf('%s minus %s',conditions{all_condition_combinations(this_cond,:)==1},conditions{all_condition_combinations(this_cond,:)==-1});
end

all_F_condition_combinations_strings = cell(1,size(all_F_condition_contrasts,1));
for this_cond = 1:size(all_F_condition_contrasts,1)
    all_F_condition_combinations_strings{this_cond} = sprintf('%s minus %s',conditions{all_F_condition_contrasts(this_cond,:)==1},conditions{all_F_condition_contrasts(this_cond,:)==-1});
end

%% Contrasts (Combined SPM for patients/controls)

contrasts = {};
cnt = 0;

if length(conditions)==2
%     cnt = cnt + 1;
%     contrasts{cnt}.name = 'Overall response (All)';
%     contrasts{cnt}.c = kron(ones(1,length(p.diagnosis_list)),[1,1,zeros(1,ncond-2)]);
%     contrasts{cnt}.type = 'F';
    
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
else %Assume standard first and no combined deviant
%     cnt = cnt + 1;
%     contrasts{cnt}.name = 'Overall response (All)';
%     contrasts{cnt}.c = kron(ones(1,length(p.diagnosis_list)),ones(1,ncond));
%     contrasts{cnt}.type = 'T';
    
    cnt = cnt + 1;
    contrasts{cnt}.name = 'Mismatch response (All)';
    contrasts{cnt}.c = kron(ones(1,length(p.diagnosis_list)),[1,repmat((-1/(ncond-1)),1,ncond-1)]);
    contrasts{cnt}.type = 'F';
    
    cnt = cnt + 1;
    contrasts{cnt}.name = ['Mismatch response controls minus allpats']; %NB weights each patient group equally regardless of size
    contrasts{cnt}.c = kron([1, -(ones(1,length(p.diagnosis_list)-1)/(length(p.diagnosis_list)-1))],[1,repmat((-1/(ncond-1)),1,ncond-1)]);
    contrasts{cnt}.type = 'T';
    
    cnt = cnt + 1;
    contrasts{cnt}.name = ['Mismatch response allpats minus controls'];
    contrasts{cnt}.c = kron([-1, (ones(1,length(p.diagnosis_list)-1)/(length(p.diagnosis_list)-1))],[1,repmat((-1/(ncond-1)),1,ncond-1)]);
    contrasts{cnt}.type = 'T';
    
    cnt = cnt + 1;
    contrasts{cnt}.name = ['Mismatch response controls vs eachgroup'];
    contrasts{cnt}.c = kron(vs_control_combinations(1:size(vs_control_combinations,1)/2,:),[1,repmat((-1/(ncond-1)),1,ncond-1)]);
    contrasts{cnt}.type = 'F';
end

% for this_cond = 1:ncond
%     cnt = cnt + 1;
%     contrasts{cnt}.name = [conditions{this_cond} '; All subjects'];
%     contrasts{cnt}.c = kron(ones(1,length(p.diagnosis_list)),[zeros(1,this_cond-1), 1, zeros(1,ncond-this_cond)]);
%     contrasts{cnt}.type = 'T';
% end

for this_grp = 1:length(p.diagnosis_list)
%     cnt = cnt + 1;
%     contrasts{cnt}.name = ['Overall response; ' p.diagnosis_list{this_grp}];
%     contrasts{cnt}.c = kron([zeros(1,this_grp-1), 1, zeros(1,length(p.diagnosis_list)-this_grp)],ones(1,ncond));
%     contrasts{cnt}.type = 'T';
    
    cnt = cnt + 1;
    contrasts{cnt}.name = ['Mismatch response; ' p.diagnosis_list{this_grp}];
    contrasts{cnt}.c = kron([zeros(1,this_grp-1), 1, zeros(1,length(p.diagnosis_list)-this_grp)],[1,repmat((-1/(ncond-1)),1,ncond-1)]);
    contrasts{cnt}.type = 'F';
    
    if this_grp~=1
    cnt = cnt + 1;
    contrasts{cnt}.name = ['Mismatch response controls minus ' p.diagnosis_list{this_grp}];
    contrasts{cnt}.c = kron([1, zeros(1,this_grp-2), -1, zeros(1,length(p.diagnosis_list)-this_grp)],[1,repmat((-1/(ncond-1)),1,ncond-1)]);
    contrasts{cnt}.type = 'T';
    
    cnt = cnt + 1;
    contrasts{cnt}.name = ['Mismatch response ' p.diagnosis_list{this_grp} ' minus controls'];
    contrasts{cnt}.c = kron([-1, zeros(1,this_grp-2), 1, zeros(1,length(p.diagnosis_list)-this_grp)],[1,repmat((-1/(ncond-1)),1,ncond-1)]);
    contrasts{cnt}.type = 'T';
    end
end

for this_cond = 1:(size(all_F_condition_contrasts,1))
    cnt = cnt + 1;
	contrasts{cnt}.name = [all_F_condition_combinations_strings{this_cond} '; All subjects F'];
    contrasts{cnt}.c = kron(ones(1,length(p.diagnosis_list)),all_F_condition_contrasts(this_cond,:));
    contrasts{cnt}.type = 'F';
end

for this_cond = 1:size(all_condition_combinations,1)
    cnt = cnt + 1;
	contrasts{cnt}.name = [all_condition_combinations_strings{this_cond} '; All subjects'];
    contrasts{cnt}.c = kron(ones(1,length(p.diagnosis_list)),all_condition_combinations(this_cond,:));
    contrasts{cnt}.type = 'T';
end


% for this_grp = 1:size(vs_control_combinations,1)
%     for this_cond = 1:size(all_condition_combinations,1)
%         cnt = cnt + 1;
%         contrasts{cnt}.name = [all_condition_combinations_strings{this_cond} '; ' vs_control_combinations_strings{this_grp}];
%         contrasts{cnt}.c = kron(vs_control_combinations(this_grp,:),all_condition_combinations(this_cond,:));
%         contrasts{cnt}.type = 'T';
%     end
% end
% 

% for this_grp = 1:size(all_group_combinations,1)
%     for this_cond = 1:size(all_condition_combinations,1)
%         cnt = cnt + 1;
%         contrasts{cnt}.name = [all_condition_combinations_strings{this_cond} '; ' all_group_combinations_strings{this_grp}];
%         contrasts{cnt}.c = kron(all_group_combinations(this_grp,:),all_condition_combinations(this_cond,:));
%         contrasts{cnt}.type = 'T';
%     end
% end
% 
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


%% Estimate models

outputfullpath = outputstem;
if ~exist(outputfullpath)
    mkdir(outputfullpath);
end

for this_subj = 1:length(Participant)
    all_groups(this_subj) = find(strcmp(Participant{this_subj}.diag,p.diagnosis_list));
end
all_group_numbers = unique(all_groups);

for this_wind = 1:length(p.time_wind_path)
    for g = 1:length(all_group_numbers)
        for this_condition = 1:length(p.these_conditions)
            this_group = all_group_numbers(g);
            these_files = squeeze(image_files(this_wind,all_groups==this_group,this_condition));
            for i = 1:length(these_files)
                files{this_group}{i}{this_condition} = these_files{i};
            end
        end
    end
    outputfullpath = [outputstem filesep p.time_wind_path{this_wind} filesep];
    if ~exist(outputfullpath)
        mkdir(outputfullpath);
    end
    
    % set up input structure for batch_spm_anova_vES
    S.imgfiles = files;
    S.maskimg = mskname;
    S.outdir = outputfullpath;
    S.contrasts = contrasts;
    S.uUFp = 1; % for M/EEG only
    
    batch_spm_anova_version_es(S); % estimate model and compute contrasts
    
end

