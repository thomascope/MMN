%%Make sure correct version of SPM running
spmpath = '/group/language/data/thomascope/spm12_fil_r6906/';
addpath(spmpath)
spm eeg

%This script ICA denoises the data for MMN analysis

script_dir = '/group/language/data/thomascope/MMN/ICA_denoise/';
%pathstem = '/imaging/tc02/Holly_MMN/ICA_denoise/';
pathstem = '/imaging/tc02/Holly_MMN/ICA_denoise_longwindow/';
data_definition_dir = '/imaging/hp02/pnfa_mmn/preprocessed/For_Thomas_dvts_sep/';
folder_structure_file_maindata = 'participant_folder_structure.m';
mridirectory = '/imaging/hp02/pnfa_mmn/preprocessed/For_Thomas_dvts_sep/mri_scans/';

addpath(script_dir)

% Define source reconstruction parameters

inv_meth        = { 'IID', 'LOR'};%,               'GS',          'LOR', 'EBB' };
time_wind_path  = { '-100_500' '150 250'};
windows         = { [-100 500], [150 250]};
wind_cnt        = 1; %Which time window for LFP extraction

p.windows = cell2mat(windows');

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

missing_scan_name = {};
missing_scan_diag = {};
missing_scan_subjnum = [];
for todonumber = 1:136
    if strcmp(Participant{todonumber}.MRI,'single_subj_T1')
        missing_scan_name{end+1} = Participant{todonumber}.name;
        missing_scan_diag{end+1} = Participant{todonumber}.diag;
        missing_scan_subjnum(end+1) = todonumber;
        if ~exist([mridirectory Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/' Participant{todonumber}.MRI '.nii'],'file')
            try
                copyfile('/group/language/data/thomascope/spm12_fil_r6906/canonical/single_subj_T1.nii',[mridirectory Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/' Participant{todonumber}.MRI '.nii'])
            end
        end
        
    end
end

% Specify parameters
%p.mod = {'MEGMAG' 'MEGPLANAR' 'EEG'}; % imaging modality (used by 'convert','convert+epoch','image','smooth','mask','firstlevel' steps) NB: EEG MUST ALWAYS BE LISTED LAST!!
p.mod = {'MEGMAG' 'MEGPLANAR'}; % imaging modality (used by 'convert','convert+epoch','image','smooth','mask','firstlevel' steps) NB: EEG MUST ALWAYS BE LISTED LAST!!
p.ref_chans = {'EOG061','EOG062','ECG063'};

% set time frequency decomposition parameters
%p.method = 'mtmconvol'; p.freqs = [30:2:90]; p.timeres = 200; p.timestep = 20; p.freqres = 30;
p.freqs = [4:2:80]; %Vector of frequencies of interest
p.method = 'morlet'; %method
p.ncycles = 7; %number of wavelet cycles
p.phase = 0; %save phase information too? (prefixed with tph)
p.tf_chans = 'All'; %cell array of channel names. Can include generic wildcards: 'All', 'EEG', 'MEG' etc.
p.timewin = [-600 1000]; %time window of interest
p.preBase_tf = -100; %TF baseline correct period with below (I don't know why this isn't a two element vector - don't blame me.)
p.postBase_tf = 0;
p.tf.method = 'LogR'; %'LogR', 'Diff', 'Rel', 'Log', 'Sqrt', 'None'
p.tf.subsample = 5; %subsample by a factor of 5 - mainly to save disk space and speed up subsequent processing. Without this, produced a file of 20-40Gb for each subject!
p.robust = 1; %Robust averaging?
p.preImageMask = -100; % pre image time (ms)
p.postImageMask = 500; % post image time (ms)
% for image smoothing
p.xSmooth = 10; % smooth for x dimension (mm)
p.ySmooth = 10; % smooth for y dimension (mm)
p.zSmooth = 10; % smooth for z (time) dimension (ms)
% define conditions
p.conditions = {'STD, DVT'};

p.contrast_labels = {'STD-DVT';'DVT-STD'};
p.contrast_weights = [1,-1;-1,1]; 

all_diagnoses = cell(1,nsubj);
for todonumber = 1:nsubj
    all_diagnoses{todonumber} = Participant{todonumber}.diag;
end
[p.diagnosis_list, ~, p.group] = unique(all_diagnoses,'stable');

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
    Participant{todonumber+old_nsubj}.namepostmerge = Participant{todonumber+old_nsubj}.name;
    Participant{todonumber+old_nsubj}.name = {};
    repeatmaxfilter = 0;
    try
        if maxfilterworkedcorrectly(todonumber) == 0
            Participant{todonumber+old_nsubj}.name = maxfilter_this_participant(thesepaths,subjfolder,this_participant_name,repeatmaxfilter)
        end
        maxfilterworkedcorrectly(todonumber) = 1
        fprintf('\n\nMaxfilter and convert complete for subject number %d,\n\n',todonumber);
    catch
        maxfilterworkedcorrectly(todonumber) = 0;
        fprintf('\n\nMaxfilter and convert failed for subject number %d\n\n',todonumber);
    end
end

%% Then copy other subjects' maxfiltered data to new folder
copycomplete = zeros(1,old_nsubj);
parfor todonumber = 1:old_nsubj
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
copyfile('/imaging/tc02/Holly_MMN/ICA_denoise/MEGArtifactTemplateTopographies.mat',[pathstem 'MEGArtifactTemplateTopographies.mat'])
copyfile('/imaging/tc02/Holly_MMN/ICA_denoise/tec_montage_all.mat',[pathstem 'tec_montage_all.mat'])
ICAcomplete = zeros(1,nsubj);
parfor todonumber = 1:nsubj
    if iscell(Participant{todonumber}.name)
        for this_file = 1:length(Participant{todonumber}.name)
            this_input_fname = [Participant{todonumber}.name{this_file} '.mat']
            this_output_folder_tail = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.namepostmerge '/']
            try
                Preprocessing_mainfunction('ICA_artifacts',this_input_fname,p,pathstem, [], this_output_folder_tail,todonumber)
                ICAcomplete(todonumber) = 1
                fprintf('\n\nICA complete for subject number %d file %d,\n\n',todonumber,this_file);
            catch
                ICAcomplete(todonumber) = 0;
                fprintf('\n\nICA failed for subject number %d file %d,\n\n',todonumber,this_file);
            end
        end
    else
        this_input_fname = [Participant{todonumber}.name '.mat']
        this_output_folder_tail = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/']
        try
            Preprocessing_mainfunction('ICA_artifacts',this_input_fname,p,pathstem, [], this_output_folder_tail,todonumber)
            ICAcomplete(todonumber) = 1
            fprintf('\n\nICA complete for subject number %d,\n\n',todonumber);
        catch
            ICAcomplete(todonumber) = 0;
            fprintf('\n\nICA failed for subject number %d,\n\n',todonumber);
        end
    end
end


%% Now run Holly's preprocessing
startagain = 0; %If want to repeat this step
Preprocesscomplete = zeros(1,nsubj);
parfor todonumber = 1:nsubj
    try
        if iscell(Participant{todonumber}.name)
            preproc_this_participant([pathstem Participant{todonumber}.groupfolder '/' Participant{todonumber}.namepostmerge '/'], Participant{todonumber}.name, startagain,p)
        else
            preproc_this_participant([pathstem Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/'], Participant{todonumber}.name, startagain,p)
        end
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
        Participant{todonumber}.name = Participant{todonumber}.namepostmerge;
    end
    megpaths = {[pathstem Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/fmraedfffM' Participant{todonumber}.name '.mat'],
                [pathstem Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/raedfffM' Participant{todonumber}.name '.mat']
                };
    mripath = [mridirectory Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/' Participant{todonumber}.MRI '.nii'];
    if ~exist(mripath,'file') && strcmp(Participant{todonumber}.MRI,'single_subj_T1')
        mripath = ['/group/language/data/thomascope/spm12_fil_r6906/canonical/single_subj_T1.nii'];
    end
    newmripath = [pathstem Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/MRI/'];
    if ~exist(newmripath,'dir')
        mkdir(newmripath)
    end
    try
        copyfile(mripath, newmripath)
    catch
        mripath = ['/group/language/data/thomascope/spm12_fil_r6906/canonical/single_subj_T1.nii'];
        warning(['Scan missing for subject ' num2str(todonumber) ' using template instead.']);
        copyfile(mripath, newmripath);
        Participant{todonumber}.MRI = 'single_subj_T1';
    end
    newmripath = [pathstem Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/MRI/' Participant{todonumber}.MRI '.nii'];
    try
        forward_model_this_subj(megpaths,newmripath, inv_meth, time_wind_path, windows)
        forwardmodelcomplete(todonumber) = 1;
        fprintf('\n\nForward modelling complete for subject number %d,\n\n',todonumber);
    catch
        forwardmodelcomplete(todonumber) = 0;
        fprintf('\n\nForward modelling failed for subject number %d\n\n',todonumber);
    end
end

%% Now extract the LFPs
LFPExtractioncomplete = zeros(1,nsubj);
parfor todonumber = 1:nsubj
    for inv_cnt = 1:length(inv_meth)
        megpaths = {[pathstem Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/' 's_' time_wind_path{wind_cnt} '_' inv_meth{inv_cnt} '_fmraedfffM' Participant{todonumber}.name '.mat'],
            [pathstem Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/' 's_' time_wind_path{wind_cnt} '_' inv_meth{inv_cnt} '_raedfffM' Participant{todonumber}.name '.mat']
            };
        outpath = [pathstem 'LFPs'];
        for thismeg = 1:length(megpaths)
            try
                
                [~] = Fullpipeline_extraction(megpaths{thismeg},Participant{todonumber}.diag,outpath,inv_cnt);
                LFPExtractioncomplete(todonumber) = LFPExtractioncomplete(todonumber)+1;
                
                
                fprintf('\n\n LFP extraction complete for subject number %d,\n\n',todonumber);
            catch
                %LFPExtractioncomplete(inv_cnt,thismeg,todonumber) = 0;
                fprintf('\n\n LFP extraction failed for subject number %d,\n\n',todonumber);
            end
        end
    end
end

%% Now plot the LFPs for sanity check
prefix = 'fmraedfffM';
val = 2; 
p.time_wind_path = time_wind_path;
p.wind_cnt = wind_cnt;
p.inv_meth = inv_meth;
p.inv_cnt = 2; %for LORETA
plot_all_LFPs(Participant,pathstem,p,prefix)

%% Now run Granger Causality and Imaginary Coherence
p.start_times = 0;
p.end_times = 500;
prefix = 'raedfffM';
decompositionworkedcorrectly = {};
for method = {'granger','coh'}
    p.method = char(method);
    decompositionworkedcorrectly{end+1} = Coherence_Connectivity(Participant,pathstem,p,prefix);
end

%% Now do a time-frequency analysis
prefix = 'raedfffM';
TFdecompositioncomplete = zeros(1,nsubj);
megpath = [];
parfor todonumber = 1:nsubj
    megpath = [pathstem Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/' prefix Participant{todonumber}.name '.mat'];
    if iscell(Participant{todonumber}.name)
        this_output_folder_tail = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.namepostmerge '/'];
    else
        this_output_folder_tail = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/'];
    end
    try
        Preprocessing_mainfunction('TF',megpath,p,pathstem, [], this_output_folder_tail,todonumber);
        TFdecompositioncomplete(todonumber) = 1;
    catch
        TFdecompositioncomplete(todonumber) = 0;
    end
end

prefix = 'tf_raedfffM';
TFaveragecomplete = zeros(1,nsubj);
megpath = [];
parfor todonumber = 1:nsubj
    megpath = [pathstem Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/' prefix Participant{todonumber}.name '.mat'];
    if iscell(Participant{todonumber}.name)
        this_output_folder_tail = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.namepostmerge '/'];
    else
        this_output_folder_tail = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/'];
    end
    try
        Preprocessing_mainfunction('average',megpath,p,pathstem, [], this_output_folder_tail,todonumber);
        TFaveragecomplete(todonumber) = 1;
    catch
        TFaveragecomplete(todonumber) = 0;
    end
end

prefix = 'mtf_raedfffM';
TFrescalecomplete = zeros(1,nsubj);
megpath = [];
parfor todonumber = 1:nsubj
    megpath = [pathstem Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/' prefix Participant{todonumber}.name '.mat'];
    if iscell(Participant{todonumber}.name)
        this_output_folder_tail = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.namepostmerge '/'];
    else
        this_output_folder_tail = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/'];
    end
    try
        Preprocessing_mainfunction('TF_rescale',megpath,p,pathstem, [], this_output_folder_tail,todonumber);
        TFrescalecomplete(todonumber) = 1;
    catch
        TFrescalecomplete(todonumber) = 0;
    end
end

prefix = 'rmtf_raedfffM';
TFweightcomplete = zeros(1,nsubj);
megpath = [];
parfor todonumber = 1:nsubj
    megpath = [pathstem Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/' prefix Participant{todonumber}.name '.mat'];
    if iscell(Participant{todonumber}.name)
        this_output_folder_tail = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.namepostmerge '/'];
    else
        this_output_folder_tail = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/'];
    end
    try
        Preprocessing_mainfunction('weight',megpath,p,pathstem, [], this_output_folder_tail,todonumber);
        TFweightcomplete(todonumber) = 1;
    catch
        TFweightcomplete(todonumber) = 0;
    end
end

prefix = 'rmtf_raedfffM';
TFimagecomplete = zeros(1,nsubj);
megpath = [];
parfor todonumber = 1:nsubj
    megpath = [pathstem Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/' prefix Participant{todonumber}.name '.mat'];
    if iscell(Participant{todonumber}.name)
        this_output_folder_tail = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.namepostmerge '/'];
    else
        this_output_folder_tail = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/'];
    end
    try
        Preprocessing_mainfunction('image',megpath,p,pathstem, [], this_output_folder_tail,todonumber);
        TFimagecomplete(todonumber) = 1;
    catch
        TFimagecomplete(todonumber) = 0;
    end
end

prefix = 'rmtf_raedfffM';
TFsmoothcomplete = zeros(1,nsubj);
megpath = [];
parfor todonumber = 1:nsubj
    megpath = [pathstem Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/' prefix Participant{todonumber}.name '.mat'];
    if iscell(Participant{todonumber}.name)
        this_output_folder_tail = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.namepostmerge '/'];
    else
        this_output_folder_tail = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/'];
    end
    try
        Preprocessing_mainfunction('smooth',megpath,p,pathstem, [], this_output_folder_tail,todonumber);
        TFsmoothcomplete(todonumber) = 1;
    catch
        TFsmoothcomplete(todonumber) = 0;
    end
end

prefix = 'rmtf_raedfffM';
TFmaskcomplete = zeros(1,1);
megpath = [];
for todonumber = 1:1 % Only need 1 mask
    megpath = [pathstem Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/' prefix Participant{todonumber}.name '.mat'];
    if iscell(Participant{todonumber}.name)
        this_output_folder_tail = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.namepostmerge '/'];
    else
        this_output_folder_tail = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/'];
    end
    try
        Preprocessing_mainfunction('mask',megpath,p,pathstem, [], this_output_folder_tail,todonumber);
        TFmaskcomplete(todonumber) = 1;
    catch
        TFmaskcomplete(todonumber) = 0;
    end
end

%% Now do TF grand averages 

prefix = 'wrmtf_raedfffM*.mat';
TFweightedgrandaveragecomplete = zeros(1,1);
for todonumber = 1:nsubj
    if iscell(Participant{todonumber}.name)
        this_output_folder_tail{todonumber}  = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.namepostmerge '/'];
    else
        this_output_folder_tail{todonumber}  = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/'];
    end
end
try
    Preprocessing_mainfunction('grand_average',prefix,p,pathstem, [], this_output_folder_tail,todonumber);
    TFweightedgrandaveragecomplete(1) = 1;
catch
    TFweightedgrandaveragecomplete(1) = 0;
end
if ~exist([pathstem 'TF_grand_averages'])
    mkdir([pathstem 'TF_grand_averages'])
end
filestomove = dir([pathstem '*weighted_grandmean*.mat']);
for i = 1:length(filestomove)
S.D = [pathstem filestomove(i).name];
S.outfile = [pathstem 'TF_grand_averages/' filestomove(i).name];
spm_eeg_copy(S)
end
delete([pathstem '*weighted_grandmean*'])

prefix = 'rmtf_raedfffM*.mat';
TFgrandaveragecomplete = zeros(1,1);
for todonumber = 1:nsubj
    if iscell(Participant{todonumber}.name)
        this_output_folder_tail{todonumber}  = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.namepostmerge '/'];
    else
        this_output_folder_tail{todonumber}  = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/'];
    end
end
try
    Preprocessing_mainfunction('grand_average',prefix,p,pathstem, [], this_output_folder_tail,todonumber);
    TFgrandaveragecomplete(1) = 1;
catch
    TFgrandaveragecomplete(1) = 0;
end
filestomove = dir([pathstem '*_grandmean*.mat']);
for i = 1:length(filestomove)
S.D = [pathstem filestomove(i).name];
S.outfile = [pathstem 'TF_grand_averages/' filestomove(i).name];
spm_eeg_copy(S)
end
delete([pathstem '*_grandmean*'])


%% Also grand average the non-TF data
prefix = 'wfmraedfffM*.mat';
weightedgrandaveragecomplete = zeros(1,1);
for todonumber = 1:nsubj
    if iscell(Participant{todonumber}.name)
        this_output_folder_tail{todonumber}  = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.namepostmerge '/'];
    else
        this_output_folder_tail{todonumber}  = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/'];
    end
end
try
    Preprocessing_mainfunction('grand_average',prefix,p,pathstem, [], this_output_folder_tail,todonumber);
    TFweightedgrandaveragecomplete(1) = 1;
catch
    TFweightedgrandaveragecomplete(1) = 0;
end

if ~exist([pathstem 'ERP_grand_averages'])
    mkdir([pathstem 'ERP_grand_averages'])
end
filestomove = dir([pathstem '*weighted_grandmean*.mat']);
for i = 1:length(filestomove)
S.D = [pathstem filestomove(i).name];
S.outfile = [pathstem 'ERP_grand_averages/' filestomove(i).name];
spm_eeg_copy(S)
end
delete([pathstem '*weighted_grandmean*'])

prefix = 'fmraedfffM*.mat';
grandaveragecomplete = zeros(1,1);
for todonumber = 1:nsubj
    if iscell(Participant{todonumber}.name)
        this_output_folder_tail{todonumber}  = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.namepostmerge '/'];
    else
        this_output_folder_tail{todonumber}  = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/'];
    end
end
try
    Preprocessing_mainfunction('grand_average',prefix,p,pathstem, [], this_output_folder_tail,todonumber);
    TFweightedgrandaveragecomplete(1) = 1;
catch
    TFweightedgrandaveragecomplete(1) = 0;
end
filestomove = dir([pathstem '*_grandmean*.mat']);
for i = 1:length(filestomove)
S.D = [pathstem filestomove(i).name];
S.outfile = [pathstem 'ERP_grand_averages/' filestomove(i).name];
spm_eeg_copy(S)
end
delete([pathstem '*_grandmean*'])


