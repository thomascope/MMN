%Make sure correct version of SPM running
spmpath = '/group/language/data/thomascope/spm12_fil_r6906/';
addpath(spmpath)
spm eeg

%This script ICA denoises the data for MMN analysis

script_dir = '/group/language/data/thomascope/MMN/ICA_denoise/';
pathstem = '/imaging/tc02/Holly_MMN/ICA_denoise/';
data_definition_dir = '/imaging/hp02/pnfa_mmn/preprocessed/For_Thomas_dvts_sep/';
folder_structure_file_maindata = 'participant_folder_structure.m';

% Define data location
run([data_definition_dir folder_structure_file_maindata]);
nsubj = size(Participant,2);

% Specify parameters
p.mod = {'MEGMAG' 'MEGPLANAR' 'EEG'}; % imaging modality (used by 'convert','convert+epoch','image','smooth','mask','firstlevel' steps) NB: EEG MUST ALWAYS BE LISTED LAST!!
p.ref_chans = {'EOG061','EOG062','ECG063'};

%Open parallel pool
try
    if ~isempty(gcp('nocreate'))
        delete(gcp)
    end
catch
    try 
        matlabpool close
    catch
    end
end

if nsubj <= 92 % Max 92 workers for 1 week
    cbupool(nsubj,'--mem-per-cpu=8G --time=167:00:00')
else
    cbupool(92,'--mem-per-cpu=8G --time=167:00:00')
end

%First copy Holly's data to new folder
todocomplete = zeros(1,nsubj);
parfor todonumber = 1:size(todoarray,1)
    this_input_full_fname = [preproc_path Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/' Participant{todonumber}.name '.mat']
    this_output_folder_tail = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/']
    try
        Preprocessing_mainfunction('Holly_data_copy',this_input_full_fname,p,pathstem, [], this_output_folder_tail,todonumber,[Participant{todonumber}.name '.mat'])
        todocomplete(todonumber) = 1
        fprintf('\n\nICA complete for subject number %d, run number %d\n\n',todoarray(todonumber,1), todoarray(todonumber,2));
    catch
        todocomplete(todonumber) = 0;
        fprintf('\n\nICA failed for subject number %d, run number %d\n\n',todoarray(todonumber,1), todoarray(todonumber,2));
    end
end

% Now run ICA_denoise
todocomplete = zeros(1,nsubj);
parfor todonumber = 1:size(todoarray,1)
    this_input_fname = [Participant{todonumber}.name '.mat']
    this_output_folder_tail = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/']
    try
        Preprocessing_mainfunction('ICA_artifacts',this_input_fname,p,pathstem, [], this_output_folder_tail,todonumber)
        todocomplete(todonumber) = 1
        fprintf('\n\nICA complete for subject number %d, run number %d\n\n',todoarray(todonumber,1), todoarray(todonumber,2));
    catch
        todocomplete(todonumber) = 0;
        fprintf('\n\nICA failed for subject number %d, run number %d\n\n',todoarray(todonumber,1), todoarray(todonumber,2));
    end
end

[preproc_path Participant{1}.groupfolder '/' Participant{1}.name '/' Participant{1}.name '.mat']
