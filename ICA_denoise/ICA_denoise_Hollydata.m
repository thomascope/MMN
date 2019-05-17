%Make sure correct version of SPM running
spmpath = '/group/language/data/thomascope/spm12_fil_r6906/';
addpath(spmpath)
spm eeg

%This script ICA denoises the data for MMN analysis

script_dir = '/group/language/data/thomascope/MMN/ICA_denoise/';
pathstem = '/imaging/tc02/Holly_MMN/ICA_denoise/';
data_definition_dir = '/imaging/hp02/pnfa_mmn/preprocessed/For_Thomas_dvts_sep/';
folder_structure_file_maindata = 'participant_folder_structure.m';

addpath(script_dir)

% Define data location
run([data_definition_dir folder_structure_file_maindata]);
nsubj = size(Participant,2);

% Specify parameters
%p.mod = {'MEGMAG' 'MEGPLANAR' 'EEG'}; % imaging modality (used by 'convert','convert+epoch','image','smooth','mask','firstlevel' steps) NB: EEG MUST ALWAYS BE LISTED LAST!!
p.mod = {'MEGMAG' 'MEGPLANAR'}; % imaging modality (used by 'convert','convert+epoch','image','smooth','mask','firstlevel' steps) NB: EEG MUST ALWAYS BE LISTED LAST!!
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
copycomplete = zeros(1,nsubj);
parfor todonumber = 1:nsubj
    this_input_full_fname = [preproc_path Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/' Participant{todonumber}.name '.mat']
    this_output_folder_tail = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/']
    try
        Preprocessing_mainfunction('Holly_data_copy',this_input_full_fname,p,pathstem, [], this_output_folder_tail,todonumber,[Participant{todonumber}.name '.mat'])
        copycomplete(todonumber) = 1
        fprintf('\n\nCopy complete for subject number %d,\n\n',todonumber);
    catch
        copycomplete(todonumber) = 0;
        fprintf('\n\nCopy failed for subject number %d\n\n',todonumber);
    end
end

% Now run ICA_denoise
ICAcomplete = zeros(1,nsubj);
parfor todonumber = 1:nsubj
    this_input_fname = [Participant{todonumber}.name '.mat']
    this_output_folder_tail = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/']
    try
        Preprocessing_mainfunction('ICA_artifacts',this_input_fname,p,pathstem, [], this_output_folder_tail,todonumber)
        ICAcomplete(todonumber) = 1
        fprintf('\n\nICA complete for subject number %d,\n\n',todonumber);
    catch
        ICAcomplete(todonumber) = 0;
        fprintf('\n\nICA failed for subject number %d\n\n',todonumber);
    end
end

% Now run Holly's preprocessing
Preprocesscomplete = zeros(1,nsubj);
parfor todonumber = 1:nsubj
    try
        preproc_this_participant([pathstem Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/'], ['M' Participant{todonumber}.name])
        Preprocesscomplete(todonumber) = 1
        fprintf('\n\nPreprocessing complete for subject number %d,\n\n',todonumber);
    catch
        Preprocesscomplete(todonumber) = 0;
        fprintf('\n\nPreprocessing failed for subject number %d\n\n',todonumber);
    end
end