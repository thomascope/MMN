%%Make sure correct version of SPM running
spmpath = '/group/language/data/thomascope/spm12_fil_r6906/';
addpath(spmpath)
spm eeg

%This script ICA denoises the data for MMN analysis

script_dir = '/group/language/data/thomascope/MMN/ICA_denoise/';
pathstem = '/imaging/tc02/Holly_MMN/ICA_denoise/';
data_definition_dir = '/imaging/hp02/pnfa_mmn/preprocessed/For_Thomas_dvts_sep/';
folder_structure_file_maindata = 'participant_folder_structure.m';
mridirectory = '/imaging/hp02/pnfa_mmn/preprocessed/For_Thomas_dvts_sep/mri_scans/';

addpath(script_dir)

% Define data location
%run([data_definition_dir folder_structure_file_maindata]);
run(folder_structure_file_maindata);

% Define MCI data location
MCI_subjs_dir = '/megdata/cbu/camcan_f/';
all_MCI_subjs = dir([MCI_subjs_dir '*pp*']);
fullpath = {};
errored_subjs = {};
MCI_fnames = {};
for i = 1:size(all_MCI_subjs)
    subdir = ls([MCI_subjs_dir all_MCI_subjs(i).name]);
    try
        fullpath{i} = ls([MCI_subjs_dir all_MCI_subjs(i).name '/' deblank(subdir) '/*mmn*.fif']);
        if ~cellfun('isempty',fullpath(i))
            MCI_fnames{end+1} = all_MCI_subjs(i).name;
        end
    catch
        errored_subjs{i} = [MCI_subjs_dir all_MCI_subjs(i).name '/' deblank(subdir) '/'];
    end
end
%All errors seem to be missing data so:
fullpath = fullpath(~cellfun('isempty',fullpath));

%Now add MCI participants to the existing path structure
%First check they aren't already there (e.g. if script run again without
%clearing workspace)
for i = 1:size(Participant,2)
    if strcmp(Participant{i}.diag,'MCI')
        error('Already some MCI patients in the Participant structure!')
    end
end
nsubj = size(Participant,2);
dicomsneedconverting = {};
nomri = {};
for todonumber = 1:size(fullpath,2)
    Participant{todonumber+nsubj}.name = MCI_fnames{todonumber};
    Participant{todonumber+nsubj}.groupfolder = 'MCI';
    Participant{todonumber+nsubj}.diag = 'MCI';
    
    try
    mri_path = dir(['/imaging/hp02/pnfa_mmn/preprocessed/For_Thomas_dvts_sep/mri_scans/MCI/' MCI_fnames{todonumber} '/*.nii']);
    
    if size(mri_path,1) ~= 1
        error(['more than one MRI found for subject ' num2str(todonumber)])
    else
        [~,mri_name,~] = fileparts(mri_path.name);
    end
    catch
        try
            mri_path = dir(['/imaging/hp02/pnfa_mmn/preprocessed/For_Thomas_dvts_sep/mri_scans/MCI/' MCI_fnames{todonumber} '/*.dcm']);
            if size(mri_path,1) ~= 1
                error(['more than one MRI found for subject ' num2str(todonumber)])
            else
                [~,mri_name,~] = fileparts(mri_path.name);
            end
            dicomsneedconverting{end+1} = mri_path;
        catch
            mri_name = 'single_subj_T1'; %for template
            nomri{end+1} = MCI_fnames{todonumber};
        end
    end
        
    Participant{todonumber+nsubj}.MRI = mri_name;
    thesepaths = strsplit(fullpath{todonumber});
    thesepaths = thesepaths(~(cellfun('isempty',thesepaths)));
    Participant{todonumber+nsubj}.MF = thesepaths;
end
disp(['Added ' num2str(size(Participant,2)-nsubj) ' MCI patients, of whom ' num2str(size(Participant,2)-nsubj-size(nomri,2)) ' have MRIs'])
old_nsubj = nsubj;
nsubj = size(Participant,2);

% Specify parameters
%p.mod = {'MEGMAG' 'MEGPLANAR' 'EEG'}; % imaging modality (used by 'convert','convert+epoch','image','smooth','mask','firstlevel' steps) NB: EEG MUST ALWAYS BE LISTED LAST!!
p.mod = {'MEGMAG' 'MEGPLANAR'}; % imaging modality (used by 'convert','convert+epoch','image','smooth','mask','firstlevel' steps) NB: EEG MUST ALWAYS BE LISTED LAST!!
p.ref_chans = {'EOG061','EOG062','ECG063'};

%% Open parallel pool
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

%% First maxfilter and convert MCI/AD data
maxfilterworkedcorrectly = zeros(1,size(fullpath,2));
assert((size(fullpath,2)+old_nsubj)==nsubj,'There don''t seem to be the right number of Participants in the definition files')
parfor todonumber = 1:size(fullpath,2)
    thesepaths = Participant{todonumber+old_nsubj}.MF
    subjfolder = [pathstem 'MCI/'];
    this_participant_name = [Participant{todonumber+old_nsubj}.name]
    try
        maxfilter_this_participant(thesepaths,subjfolder,this_participant_name)
        maxfilterworkedcorrectly(todonumber) = 1
        fprintf('\n\nCopy complete for subject number %d,\n\n',todonumber);
    catch
        maxfilterworkedcorrectly(todonumber) = 0;
        fprintf('\n\nCopy failed for subject number %d\n\n',todonumber);
    end
end

%% Then copy other subjects' maxfiltered data to new folder
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

%% Now run ICA_denoise
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

%% Now run Holly's preprocessing
startagain = 1; %If want to repeat this step
Preprocesscomplete = zeros(1,nsubj);
parfor todonumber = 1:nsubj
    try
        preproc_this_participant([pathstem Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/'], ['M' Participant{todonumber}.name],startagain)
        Preprocesscomplete(todonumber) = 1
        fprintf('\n\nPreprocessing complete for subject number %d,\n\n',todonumber);
    catch
        Preprocesscomplete(todonumber) = 0;
        fprintf('\n\nPreprocessing failed for subject number %d\n\n',todonumber);
    end
end

%% Now specify forward model
forwardmodelcomplete = zeros(1,nsubj);
parfor todonumber = 1:nsubj
    try
        forward_model_this_subj({[pathstem Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/fmraedfffM' Participant{todonumber}.name '.mat']},{[mridirectory Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/' Participant{todonumber}.MRI '.nii']})
        forwardmodelcomplete(todonumber) = 1
        fprintf('\n\Forward modelling complete for subject number %d,\n\n',todonumber);
    catch
        forwardmodelcomplete(todonumber) = 0;
        fprintf('\n\nForward modelling failed for subject number %d\n\n',todonumber);
    end
end
forward_model_this_subj(megpath, mripath)
