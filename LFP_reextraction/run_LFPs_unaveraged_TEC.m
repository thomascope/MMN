%% NB: THIS DOES NOT WORK AS THERE IS NO RAW FILE - CLEANED UP
n
clear all
clc
%addpath(genpath('/imaging/rowe/archive/users/hp02/spm12b'));
addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r6906'));
addpath(genpath('/imaging/local/software/mne'));
addpath(genpath('/imaging/rowe/archive/users/hp02/mmn_08/analysis_spm/new_spm_functions'));
% Root directory for EMEG data
bwd = '';

addpath('/imaging/rowe/archive/users/hp02/pnfa_mmn/maxfilter')
% subject paths to names:
matched_controls_fif_paths;
patient_filenames;
subjects = [match_HCs_names; ftd_folder; pca_folder; vespa_folder];

sourceloc_path = '/imaging/mlr/users/tc02/scratch/LFPs_TA_unaveraged';
if ~exist(sourceloc_path,'dir')
    mkdir(sourceloc_path)
end
preproc_path1 = '/imaging/rowe/archive/users/hp02/pnfa_mmn/preprocessed/For_Thomas_dvts_sep/matched_HCs';
preproc_path2 = '/imaging/rowe/archive/users/hp02/pnfa_mmn/preprocessed/For_Thomas_dvts_sep';
%% open matlabpool if required
%ParType = 0;  % Fun on Login machines (not generally advised!)
%ParType = 1;   % Run maxfilter call on Compute machines using spmd (faster)
% ParType = 2;   % Run on multiple Compute machines using parfar (best, but less feedback if crashes)
%
%
% if ParType
%     if matlabpool('size')==0;
%         %MaxNsubs = 1;
%         %if ParType == 2
%         %    for g=1:length(cbu_codes)
%         %        MaxNsubs = max([MaxNsubs length(cbu_codes{g})]);
%         %    end
%         %end
%         P = cbupool(96);
%         matlabpool(P);
%     end
% end


%% Define subjs' preprocessed files

prefix = 'raedfff';


for i = 1:length(match_HCs)
    subjs{i} = sprintf('%s%s/%s%s.mat',preproc_path1, match_HCs_names{i}, prefix, match_HCs_names{i});   % letter strings for subject-specific paths
    subj_dotdat{i} = sprintf('%s%s/%s%s.dat',preproc_path1, match_HCs_names{i}, prefix, match_HCs_names{i});
    MRI_dir{i} = strcat('/imaging/rowe/archive/users/hp02/pnfa_mmn/mri_scans/matched_HCs/', match_HCs_names{i}, '/');  % directory with this participant's MRI (*.img)
    subj_folder{i} = 'matched_HCs';
end

for i = length(match_HCs)+1:length(subjects)
    subjs{i} = sprintf('%s/%s',preproc_path2, all_preproc_path_mat{i-length(match_HCs)});   % letter strings for subject-specific paths
    subj_dotdat{i} = sprintf('%s/',preproc_path2, all_preproc_path_dat{i-length(match_HCs)});
    MRI_dir{i} = strcat('/imaging/rowe/archive/users/hp02/pnfa_mmn/mri_scans/', all_mri_paths{i-length(match_HCs)});  % directory with this participant's MRI (*.img)
    MRI_ydir{i} = strcat('/imaging/rowe/archive/users/hp02/pnfa_mmn/mri_scans/', all_mri_ypaths{i-length(match_HCs)});  % directory with this participant's MRI (*.img)
    
end
for i = length(match_HCs)+ 1:length(match_HCs)+23; subj_folder{i} = 'bvFTD'; end
for i = length(subj_folder)+ 1:length(subj_folder)+15; subj_folder{i} = 'pca'; end
for i = length(subj_folder)+ 1:length(subj_folder)+10; subj_folder{i} = 'pnfa'; end

nr_subs = length(subjects);
fprintf(1, 'Going to process %d data sets\n', nr_subs);

%% Paths to forward modelled data

prefix = 'raedfff';
for_path = '/imaging/rowe/archive/users/hp02/pnfa_mmn/forward_modelling/MSP';

for i = 1:length(match_HCs)
    subjs_for{i} = sprintf('%s/matched_HCs/%s/%s%s.mat',for_path, match_HCs_names{i}, prefix, match_HCs_names{i});   % letter strings for subject-specific paths
    subj_for_dotdat{i} = sprintf('%s/matched_HCs/%s/%s%s.dat',for_path, match_HCs_names{i}, prefix, match_HCs_names{i});
    
end

for i = length(match_HCs)+1:length(subjects)
    subjs_for{i} = sprintf('%s/%s',for_path, all_preproc_path_mat{i-length(match_HCs)});   % letter strings for subject-specific paths
    subj_for_dotdat{i} = sprintf('%s/%s',for_path, all_preproc_path_dat{i-length(match_HCs)});
    
end

Inversionworkedcorrectly = zeros(1,length(subjs_for));
Extractionworkedcorrectly = zeros(1, length(subjs_for));
%% Inversion and LFP extraction

parfor ss = 1:length(subjs_for)
    try
        ss
        %clear D folder
        D = subjs_for{ss};
        folder = subj_folder{ss};
        
        [~,D] = DoInvTA_TEC(D,folder,sourceloc_path,2); %set this to be the 2nd inversion - LORETA
        Inversionworkedcorrectly(ss) = 1;
        try
            D = DoExtractmTA_TEC(D,folder,sourceloc_path,2);
            %             S.D = [D.path '/' D.fname];
            %             S.outfile = ['/imaging/mlr/users/tc02/scratch/LFPs_TA/All_Ds/D_' num2str(ss) '.mat'];
            %             spm_eeg_copy(S)
            Extractionworkedcorrectly(ss) = 1;
        catch
            try
                [~,D] = DoInvTA_noEEG_TEC(D,folder,sourceloc_path,2); %set this to be the 2nd inversion - LORETA
                Inversionworkedcorrectly(ss) = 1;
                try
                    D = DoExtractmTA_TEC(D,folder,sourceloc_path,2);
                    %             S.D = [D.path '/' D.fname];
                    %             S.outfile = ['/imaging/mlr/users/tc02/scratch/LFPs_TA/All_Ds/D_' num2str(ss) '.mat'];
                    %             spm_eeg_copy(S)
                    Extractionworkedcorrectly(ss) = 1;
                catch
                    Extractionworkedcorrectly(ss) = 0;
                end
            catch
                Inversionworkedcorrectly(ss) = 0;
            end
        end
    catch
        Inversionworkedcorrectly(ss) = 0;
    end
    
end
    %50 failed
    
    % if ParType
    %     matlabpool close force CBU_Cluster
    % end
