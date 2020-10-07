%%Make sure correct version of SPM running
spmpath = '/group/language/data/thomascope/spm12_fil_r6906/';
thisspm = which('spm');
if ~strcmp(thisspm(1:end-5), spmpath)
    rmpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/'));
    addpath(spmpath)
    spm eeg
end

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

p.time_wind_path = time_wind_path;
p.inv_meth = inv_meth;
p.windows = cell2mat(windows');

% Define data location
%run([data_definition_dir folder_structure_file_maindata]);
run(folder_structure_file_maindata);

% Define MCI data location
biomarker_positive_mci = {'meg17_0144_pp110119'
'meg18_0065_pp113409'
'meg17_0199_pp114823'
'meg17_0217_pp117582'
'meg17_0171_pp119159'
'meg17_0196_pp126264'
'meg18_0011_pp167931'
'meg18_0066_pp167967'
'meg17_0193_pp175738'
'meg17_0163_pp183667'
'meg18_0001_pp187628'
'meg17_0153_pp196609'
'meg17_0116_pp142632'
'meg17_0159_pp108210'
'meg17_0247_pp114097'
};

biomarker_negative_mci = {'meg17_0200_pp111738'
'meg18_0051_pp112035'
'meg17_0216_pp128346'
'meg17_0162_pp135832'
'meg17_0206_pp136246'
'meg17_0160_pp167487'
'meg18_0053_pp137551'
};

biomarker_unknown_mci = {'meg17_0154_pp136072'
'meg17_0248_pp156841'
'meg17_0108_pp167844'
'meg17_0120_pp168080'
'meg18_0047_pp170827'
'meg17_0117_pp105571'
};

biomarker_unknown_AD = {'meg18_0010_pp102319'
'meg18_0039_pp113615'
'meg17_0238_pp117411'
'meg18_0048_pp138368'
'meg17_0240_pp141038'
'meg18_0052_pp142409'
'meg18_0035'
'meg18_0003_pp153538'
'meg18_0040_pp155559'
'meg17_0246_pp176327'
'meg18_0042_pp183367'
'meg17_0136_pp185442'
'meg18_0041_pp196451'
};

MCI_subjs_dir = '/megdata/cbu/camcan_f/';
all_MCI_subjs = dir([MCI_subjs_dir '*pp*']);
fullpath = {};
errored_subjs = {};
MCI_fnames = {};
for i = 1:size(all_MCI_subjs)
%     if any(strcmp(all_MCI_subjs(i).name,biomarker_negative_mci)) % exclude biomarker negative
%         continue
%     end
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
biomarker_positive_mci_indices = [];
biomarker_negative_mci_indices = [];
biomarker_unknown_mci_indices = [];
biomarker_unknown_AD_indices = [];
for todonumber = 1:size(fullpath,2)
    Participant{todonumber+nsubj}.name = MCI_fnames{todonumber};
    Participant{todonumber+nsubj}.groupfolder = 'MCI';
    Participant{todonumber+nsubj}.diag = 'MCI';
    
    if any(strcmp(MCI_fnames{todonumber},biomarker_positive_mci))
        biomarker_positive_mci_indices(end+1) = todonumber+nsubj;
    end
    if any(strcmp(MCI_fnames{todonumber},biomarker_negative_mci))
        biomarker_negative_mci_indices(end+1) = todonumber+nsubj;
    end
    if any(strcmp(MCI_fnames{todonumber},biomarker_unknown_mci))
        biomarker_unknown_mci_indices(end+1) = todonumber+nsubj;
    end
    if any(strcmp(MCI_fnames{todonumber},biomarker_unknown_AD))
        biomarker_unknown_AD_indices(end+1) = todonumber+nsubj;
    end

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
%p.conditions = {'STD, DVT'};
%p.conditions = {'STD','DVT','Loc','Int','Dur','Gap','Freq','Loc_L','Freq_hi','Int_hi','Loc_R','Freq_lo','Int_lo'};
p.conditions = {'STD','DVT','location','intensity','duration','gap','frequency','location_L','frequency_high','intensity_high','location_R','frequency_low','intensity_low'};
p.preBase = -100; %TF baseline correct period with below (I don't know why this isn't a two element vector - don't blame me.)
p.postBase = 0;

p.contrast_labels = {'STD-DVT';'DVT-STD';'STD-Loc';'STD-Int';'STD-Dur';'STD-Gap';'STD-Freq'};
p.contrast_weights = [1,-1,0,0,0,0,0;-1,1,0,0,0,0,0;1,0,-1,0,0,0,0;1,0,0,-1,0,0,0;1,0,0,0,-1,0,0;1,0,0,0,0,-1,0;1,0,0,0,0,0,-1]; 

%% Exclude MCIs with unknown or negative biomarker status
for i = biomarker_positive_mci_indices
    Participant{i}.diag = 'ADMCI';
    Participant{i}.groupfolder = 'ADMCI';
end
for i = biomarker_negative_mci_indices
    Participant{i}.diag = 'MCI_neg';
end
for i = biomarker_unknown_mci_indices
    Participant{i}.diag = 'MCI_unk';
end
for i = biomarker_unknown_AD_indices
    Participant{i}.diag = 'ADMCI';
    Participant{i}.groupfolder = 'ADMCI';
end

Participant([biomarker_negative_mci_indices, biomarker_unknown_mci_indices]) = [];
nsubj = nsubj - (length(biomarker_negative_mci_indices)+length(biomarker_unknown_mci_indices));

all_diagnoses = cell(1,nsubj);
all_names = cell(1,nsubj);
for todonumber = 1:nsubj
    all_diagnoses{todonumber} = Participant{todonumber}.diag;
    all_names{todonumber} = Participant{todonumber}.name;
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
    cbupool(nsubj,'--mem-per-cpu=8G --time=167:00:00 --exclude=node-i[01-15]')
else
    cbupool(92,'--mem-per-cpu=8G --time=167:00:00 --exclude=node-i[01-15]')
end

% %% First maxfilter and convert MCI/AD data
% maxfilterworkedcorrectly = zeros(1,size(fullpath,2));
% assert((size(fullpath,2)+old_nsubj)==nsubj,'There don''t seem to be the right number of Participants in the definition files')
% parfor todonumber = 1:size(fullpath,2)
%     thesepaths = Participant{todonumber+old_nsubj}.MF
%     subjfolder = [pathstem 'MCI/'];
%     this_participant_name = [Participant{todonumber+old_nsubj}.name]
%     Participant{todonumber+old_nsubj}.namepostmerge = Participant{todonumber+old_nsubj}.name;
%     Participant{todonumber+old_nsubj}.name = {};
%     repeatmaxfilter = 0;
%     try
%         if maxfilterworkedcorrectly(todonumber) == 0
%             Participant{todonumber+old_nsubj}.name = maxfilter_this_participant(thesepaths,subjfolder,this_participant_name,repeatmaxfilter)
%         end
%         maxfilterworkedcorrectly(todonumber) = 1
%         fprintf('\n\nMaxfilter and convert complete for subject number %d,\n\n',todonumber);
%     catch
%         maxfilterworkedcorrectly(todonumber) = 0;
%         fprintf('\n\nMaxfilter and convert failed for subject number %d\n\n',todonumber);
%     end
% end
% 
% %% Then copy other subjects' maxfiltered data to new folder
% copycomplete = zeros(1,old_nsubj);
% parfor todonumber = 1:old_nsubj
%     this_input_full_fname = [preproc_path Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/' Participant{todonumber}.name '.mat']
%     this_output_folder_tail = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/']
%     try
%         Preprocessing_mainfunction('Holly_data_copy',this_input_full_fname,p,pathstem, [], this_output_folder_tail,todonumber,[Participant{todonumber}.name '.mat'])
%         copycomplete(todonumber) = 1
%         fprintf('\n\nCopy complete for subject number %d,\n\n',todonumber);
%     catch
%         copycomplete(todonumber) = 0;
%         fprintf('\n\nCopy failed for subject number %d\n\n',todonumber);
%     end
% end

%% ensure naming convention is as desired, if maxfilter step has not run
for todonumber = 1:nsubj
    if iscell(Participant{todonumber}.MF) && size(Participant{todonumber}.MF,2) > 1 && ~iscell(Participant{todonumber}.name)
        this_participant_name = [Participant{todonumber}.name];
        Participant{todonumber}.namepostmerge = Participant{todonumber}.name;
        Participant{todonumber}.name = {};
        for i = 1:size(Participant{todonumber}.MF,2)
            Participant{todonumber}.name{i} = [this_participant_name '_' num2str(i)];
        end
    end
end

%% Copy maxfiltered data to new directory structure
preproc_path_tc = '/imaging/tc02/Holly_MMN/ICA_denoise/';
copycomplete = zeros(1,old_nsubj);
parfor todonumber = 1:nsubj
    if iscell(Participant{todonumber}.name)
        for this_file = 1:length(Participant{todonumber}.name)
            this_input_full_fname = [preproc_path_tc Participant{todonumber}.groupfolder '/' Participant{todonumber}.namepostmerge '/' Participant{todonumber}.name{this_file} '.mat']
            this_output_folder_tail = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.namepostmerge '/']
            try
                Preprocessing_mainfunction('Holly_data_copy',this_input_full_fname,p,pathstem, [], this_output_folder_tail,todonumber,[Participant{todonumber}.name{this_file} '.mat'])
                copycomplete(todonumber) = 1
                fprintf('\n\nCopy complete for subject number %d,\n\n',todonumber);
            catch
                copycomplete(todonumber) = 0;
                fprintf('\n\nCopy failed for subject number %d\n\n',todonumber);
            end
        end
    else
        
        this_input_full_fname = [preproc_path_tc Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/' Participant{todonumber}.name '.mat']
        this_output_folder_tail = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/']
        try
            %Preprocessing_mainfunction('Holly_data_copy',this_input_full_fname,p,pathstem, [], this_output_folder_tail,todonumber,[Participant{todonumber}.name '.mat'])
            copycomplete(todonumber) = 1
            fprintf('\n\nCopy complete for subject number %d,\n\n',todonumber);
        catch
            copycomplete(todonumber) = 0;
            fprintf('\n\nCopy failed for subject number %d\n\n',todonumber);
        end
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
startagain = 1; %If want to repeat this step
Preprocesscomplete = zeros(1,nsubj);
parfor todonumber = 1:nsubj
    try
        if Preprocesscomplete(todonumber)~=1
        if iscell(Participant{todonumber}.name)
            preproc_this_participant([pathstem Participant{todonumber}.groupfolder '/' Participant{todonumber}.namepostmerge '/'], Participant{todonumber}.name, startagain,p)
        else
            preproc_this_participant([pathstem Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/'], Participant{todonumber}.name, startagain,p)
        end
        Preprocesscomplete(todonumber) = 1
        fprintf('\n\nPreprocessing complete for subject number %d,\n\n',todonumber);
        end
    catch
        Preprocesscomplete(todonumber) = 0;
        fprintf('\n\nPreprocessing failed for subject number %d\n\n',todonumber);
    end
end


%% Now specify forward model
forwardmodelcomplete = zeros(1,nsubj);
parfor todonumber = 1:nsubj
    if forwardmodelcomplete(todonumber)~=1
        
        try
            Participant{todonumber}.name = Participant{todonumber}.namepostmerge;
        end
        megpaths = {[pathstem Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/fmbraedfffM' Participant{todonumber}.name '.mat'],
            [pathstem Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/braedfffM' Participant{todonumber}.name '.mat']
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
            forward_model_this_subj(megpaths,newmripath, p)
            forwardmodelcomplete(todonumber) = 1;
            fprintf('\n\nForward modelling complete for subject number %d,\n\n',todonumber);
        catch
            forwardmodelcomplete(todonumber) = 0;
            fprintf('\n\nForward modelling failed for subject number %d\n\n',todonumber);
        end
    end
end

%% Now extract the LFPs
LFPExtractioncomplete = zeros(1,nsubj);
parfor todonumber = 1:nsubj
    if LFPExtractioncomplete(todonumber)~=4
        for inv_cnt = 1:length(inv_meth)
            megpaths = {[pathstem Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/' 's_' time_wind_path{wind_cnt} '_' inv_meth{inv_cnt} '_fmbraedfffM' Participant{todonumber}.name '.mat'],
                [pathstem Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/' 's_' time_wind_path{wind_cnt} '_' inv_meth{inv_cnt} '_braedfffM' Participant{todonumber}.name '.mat']
                };
            outpath = [pathstem 'LFPs'];
            for thismeg = 1:length(megpaths)
                try
                    
                    [~] = Fullpipeline_extraction(megpaths{thismeg},Participant{todonumber}.diag,outpath,time_windowt);
                    LFPExtractioncomplete(todonumber) = LFPExtractioncomplete(todonumber)+1;
                    
                    
                    fprintf('\n\n LFP extraction complete for subject number %d,\n\n',todonumber);
                catch
                    %LFPExtractioncomplete(inv_cnt,thismeg,todonumber) = 0;
                    fprintf('\n\n LFP extraction failed for subject number %d,\n\n',todonumber);
                end
            end
        end
    end
end

%% Now baseline correct the LFPs
LFPBaselinecomplete = zeros(1,nsubj);
parfor todonumber = 1:nsubj
    if LFPBaselinecomplete(todonumber)~=4
        for inv_cnt = 1:length(inv_meth)
            this_input_fname = {['8LFP_s_' time_wind_path{wind_cnt} '_' inv_meth{inv_cnt} '_fmbraedfffM' Participant{todonumber}.name '.mat'],
                ['8LFP_s_' time_wind_path{wind_cnt} '_' inv_meth{inv_cnt} '_braedfffM' Participant{todonumber}.name '.mat']
                };
            this_output_folder_tail = [Participant{todonumber}.diag '/']
            for thismeg = 1:length(this_input_fname)
                try
                    Preprocessing_mainfunction('baseline',this_input_fname{thismeg},p,[pathstem 'LFPs/'], [], this_output_folder_tail,todonumber)
                    LFPBaselinecomplete(todonumber) = LFPBaselinecomplete(todonumber) + 1;
                    if LFPBaselinecomplete(todonumber) == 4
                        fprintf('\n\nLFP Baseline complete for subject number %d,\n\n',todonumber);
                    end
                catch
                    LFPBaselinecomplete(todonumber) = 0;
                    fprintf('\n\nLFP Baseline failed for subject number %d,\n\n',todonumber);
                end
            end
        end
    end
end

%% Now plot the LFPs for sanity check
for todonumber = 1:nsubj
        try
        Participant{todonumber}.name = Participant{todonumber}.namepostmerge;
        end
end
prefix = 'fmbraedfffM';
val = 2; %for LORETA
%val = 1 %for IID
p.time_wind_path = time_wind_path;
p.wind_cnt = wind_cnt;
p.inv_meth = inv_meth;
p.inv_cnt = val;
baselined = 1;
%plot_all_LFPs(Participant,pathstem,p,prefix)
plot_MMN_bytype_LFP(Participant,pathstem,p,prefix,baselined)
% 
% %% Now plot the baseline corrected LFPs for sanity check
% prefix = 'bfmbraedfffM';
% val = 2; 
% p.time_wind_path = time_wind_path;
% p.wind_cnt = wind_cnt;
% p.inv_meth = inv_meth;
% p.inv_cnt = 2; %for LORETA
% plot_all_LFPs(Participant,pathstem,p,prefix)

%% Now run Granger Causality and Imaginary Coherence
p.start_times = 0;
p.end_times = 500;
prefix = 'braedfffM';
decompositionworkedcorrectly = {};
val = 2; %for LORETA
%val = 1 %for IID
p.time_wind_path = time_wind_path;
p.wind_cnt = wind_cnt;
p.inv_meth = inv_meth;
p.inv_cnt = val; 

for todonumber = 1:nsubj
        try
        Participant{todonumber}.name = Participant{todonumber}.namepostmerge;
        end
end
for method = {'granger','coh'}
    p.decompmethod = char(method);
    decompositionworkedcorrectly{end+1} = Coherence_Connectivity(Participant,pathstem,p,prefix);
end

%% Now do a time-frequency analysis
for todonumber = 1:nsubj
        try
        Participant{todonumber}.name = Participant{todonumber}.namepostmerge;
        end
end

prefix = 'braedfffM';
TFdecompositioncomplete = zeros(1,nsubj);
megpath = [];
parfor todonumber = 1:nsubj
    %if todonumber ~= 101 && todonumber ~= 102 && todonumber ~= 135 && todonumber ~= 136 ; continue; end
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

prefix = 'tf_braedfffM';
TFaveragecomplete = zeros(1,nsubj);
megpath = [];
parfor todonumber = 1:nsubj
        %if todonumber ~= 101 && todonumber ~= 102 && todonumber ~= 135 ; continue; end

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

prefix = 'mtf_braedfffM';
TFrescalecomplete = zeros(1,nsubj);
megpath = [];
parfor todonumber = 1:nsubj
        %if todonumber ~= 101 && todonumber ~= 102 && todonumber ~= 135 ; continue; end

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

prefix = 'rmtf_braedfffM';
TFweightcomplete = zeros(1,nsubj);
megpath = [];
parfor todonumber = 1:nsubj
        %if todonumber ~= 101 && todonumber ~= 102 && todonumber ~= 135 ; continue; end

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

prefix = 'rmtf_braedfffM';
TFimagecomplete = zeros(1,nsubj);
megpath = [];
parfor todonumber = 1:nsubj
        %if todonumber ~= 101 && todonumber ~= 102 && todonumber ~= 135 ; continue; end

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

prefix = 'rmtf_braedfffM';
TFsmoothcomplete = zeros(1,nsubj);
megpath = [];
parfor todonumber = 1:nsubj
        %if todonumber ~= 101 && todonumber ~= 102 && todonumber ~= 135 ; continue; end

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

prefix = 'rmtf_braedfffM';
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

prefix = 'wrmtf_braedfffM*.mat';e
TFweightedgrandaveragecomplete = zeros(1,1);
this_output_folder_tail = {};
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

prefix = 'rmtf_braedfffM*.mat';
TFgrandaveragecomplete = zeros(1,1);
this_output_folder_tail = {};
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

%% Now do second level analysis on the TF data
prefix = 'rmtf_braedfffM';
TFsecondlevelcomplete = zeros(1,1);
this_output_folder_tail = {};
p.mod = {'MEGMAG', 'MEGPLANAR'};
for todonumber = 1:nsubj
    if iscell(Participant{todonumber}.name)
        this_output_folder_tail{todonumber}  = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.namepostmerge '/'];
    else
        this_output_folder_tail{todonumber}  = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/'];
    end
end
p.all_conditions = p.conditions;
p.conditions = {'STD','DVT'};
try
    tc_batch_SPM(prefix,this_output_folder_tail,pathstem,p);
    secondlevelcomplete(1) = 1;
catch
    secondlevelcomplete(1) = 0;
end
p.conditions = {'STD','location','intensity','duration','gap','frequency'};
try
    tc_batch_SPM(prefix,this_output_folder_tail,pathstem,p);
    secondlevelcomplete(1) = 1;
catch
    secondlevelcomplete(1) = 0;
end
p.conditions = p.all_conditions;


%% Also grand average the non-TF data
prefix = 'wfmbraedfffM*.mat';
weightedgrandaveragecomplete = zeros(1,1);
this_output_folder_tail = {};
for todonumber = 1:nsubj
    if iscell(Participant{todonumber}.name)
        this_output_folder_tail{todonumber}  = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.namepostmerge '/'];
    else
        this_output_folder_tail{todonumber}  = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/'];
    end
end
try
    Preprocessing_mainfunction('grand_average',prefix,p,pathstem, [], this_output_folder_tail,todonumber);
    weightedgrandaveragecomplete(1) = 1;
catch
    weightedgrandaveragecomplete(1) = 0;
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

prefix = 'fmbraedfffM*.mat';
grandaveragecomplete = zeros(1,1);
this_output_folder_tail = {};
for todonumber = 1:nsubj
    if iscell(Participant{todonumber}.name)
        this_output_folder_tail{todonumber}  = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.namepostmerge '/'];
    else
        this_output_folder_tail{todonumber}  = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/'];
    end
end
try
    Preprocessing_mainfunction('grand_average',prefix,p,pathstem, [], this_output_folder_tail,todonumber);
    grandaveragecomplete(1) = 1;
catch
    grandaveragecomplete(1) = 0;
end
filestomove = dir([pathstem '*_grandmean*.mat']);
for i = 1:length(filestomove)
S.D = [pathstem filestomove(i).name];
S.outfile = [pathstem 'ERP_grand_averages/' filestomove(i).name];
spm_eeg_copy(S)
end
delete([pathstem '*_grandmean*'])

%% Create images for statistical analysis on the non-TF data

p.mod = {'MEGMAG', 'MEGCOMB'};
prefix = 'PfmbraedfffM';
imagecomplete = zeros(1,nsubj);
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
        imagecomplete(todonumber) = 1;
    catch
        imagecomplete(todonumber) = 0;
    end
end

prefix = 'PfmbraedfffM';
smoothcomplete = zeros(1,nsubj);
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
        smoothcomplete(todonumber) = 1;
    catch
        smoothcomplete(todonumber) = 0;
    end
end

prefix = 'PfmbraedfffM';
maskcomplete = zeros(1,1);
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
        maskcomplete(todonumber) = 1;
    catch
        maskcomplete(todonumber) = 0;
    end
end

%% Now do second level analysis on the non-TF data

prefix = 'PfmbraedfffM';
secondlevelcomplete = zeros(1,1);
this_output_folder_tail = {};
p.mod = {'MEGMAG', 'MEGCOMB'};
for todonumber = 1:nsubj
    if iscell(Participant{todonumber}.name)
        this_output_folder_tail{todonumber}  = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.namepostmerge '/'];
    else
        this_output_folder_tail{todonumber}  = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/'];
    end
end
p.all_conditions = p.conditions;
p.conditions = {'STD','DVT'};
try
    tc_batch_SPM(prefix,this_output_folder_tail,pathstem,p);
    secondlevelcomplete(1) = 1;
catch
    secondlevelcomplete(1) = 0;
end
p.conditions = {'STD','location','intensity','duration','gap','frequency'};
try
    tc_batch_SPM(prefix,this_output_folder_tail,pathstem,p);
    secondlevelcomplete(1) = 1;
catch
    secondlevelcomplete(1) = 0;
end
p.conditions = p.all_conditions;

%% Now run Tallie's extended DCM
p.start_times = 0;
p.end_times = 400;
prefix = 'fmbraedfffM';
val = 2; %for LORETA
%val = 1 %for IID
p.time_wind_path = time_wind_path;
p.wind_cnt = wind_cnt;
p.inv_meth = inv_meth;
p.inv_cnt = val;
conditions_to_invert = {'STD','DVT','location','intensity','duration','gap','frequency'};
%conditions_to_invert = {'STD','DVT'}; % Just do standards and deviants for now
%conditions_to_invert = {'location','intensity','duration','gap','frequency'};

%Open a parallel pool with lots of memory and spmd disabled to allow
%continuation if a worker fails
Poolinfo = cbupool(48,'--mem-per-cpu=16G --time=167:00:00 --exclude=node-i[01-15]');
parpool(Poolinfo,Poolinfo.NumWorkers,'SpmdEnabled',false);

clear all_names
for i = 1:length(Participant)
all_names{i} = Participant{i}.name;
end

% Parallelise subject and condition to avoid failure stoppages
all_condition_numbers = 1:length(conditions_to_invert);
all_subjects = 1:nsubj;
allrunsarray = [];
allrunsarray=combvec(all_subjects,all_condition_numbers)';
extDCMcomplete = zeros(1,size(allrunsarray,1));
p.subjcntforcondition = 1;
p.conditions = conditions_to_invert;
p.multilevel = 0; %for first run
parfor todonumber = 1:size(allrunsarray,1)
    this_input_fname = {['b8LFP_s_' time_wind_path{wind_cnt} '_' inv_meth{p.inv_cnt} '_' prefix Participant{allrunsarray(todonumber,1)}.name '.mat']};
    this_output_folder_tail = [Participant{allrunsarray(todonumber,1)}.diag '/']
    pause(mod(todonumber,60)); %Introduce a pause to stagger the workers - otherwise sometimes the pool fails if trying to read or write simultaneously 
    for thismeg = 1:length(this_input_fname)
        try
            Preprocessing_mainfunction('extDCM',this_input_fname{thismeg},p,[pathstem 'LFPs/'], [], this_output_folder_tail,allrunsarray(todonumber,2))
            extDCMcomplete(todonumber) = extDCMcomplete(todonumber) + 1;
            fprintf('\n\nLFP DCM modelling complete for run number %d,\n\n',todonumber);
        catch
            extDCMcomplete(todonumber) = 0;
            fprintf('\n\nLFP DCM modelling failed for run number %d,\n\n',todonumber);
        end
    end
end

% Now repeat for those few subjects who failed integration, using the posterior as a prior
extDCM_directory = '/imaging/tc02/Holly_MMN/extDCMs/';
%conditions_to_invert = {'STD','DVT','location','intensity','duration','gap','frequency'};
p.conditions = conditions_to_invert;
[subjcondpair] = find_failed_extDCM_integrations(extDCM_directory,conditions_to_invert,Participant);

these_STDs = {}; %Now find those subjects where the standard failed to integrate
for i = 1:size(subjcondpair,1)
if strcmp(subjcondpair{i,2},'STD')
these_STDs{end+1} = subjcondpair{i,1};
end
end

these_DVTs = {};
these_others = [];
for i = 1:size(subjcondpair,1)
    if strcmp(subjcondpair{i,2},'DVT') && any(contains(these_STDs,subjcondpair{i,1}))
        these_DVTs{end+1} = subjcondpair{i,1};
    elseif ~strcmp(subjcondpair{i,2},'STD') && any(contains(these_STDs,subjcondpair{i,1})) 
        these_others(end+1) = i;
    end
end

if ~isempty(these_DVTs)
    error('Both the STD and the DVT integration contain NaNs for at least one subject. This is fatally problematic for loading the posteriors as priors')
end

firstwavedata = subjcondpair(setdiff(1:size(subjcondpair,1),these_others),:);
secondwavedata = subjcondpair(intersect(1:size(subjcondpair,1),these_others),:);

subjcondpair = firstwavedata;
p.subjcntforcondition = 1;
p.multilevel = 1;
parfor this_one = 1:size(subjcondpair,1)
    this_input_fname = {['b8LFP_s_' time_wind_path{wind_cnt} '_' inv_meth{p.inv_cnt} '_' prefix subjcondpair{this_one,1} '.mat']};
    this_output_folder_tail = [Participant{find(strcmp(subjcondpair{this_one,1},all_names))}.diag '/']
    this_cond = find(strcmp(subjcondpair{this_one,2}, p.conditions));
    %pause(mod(this_one*30,90)); %Introduce a pause to stagger the workers - otherwise sometimes the pool fails if trying to read or write simultaneously 
    for thismeg = 1:length(this_input_fname)
        try
            Preprocessing_mainfunction('extDCM',this_input_fname{thismeg},p,[pathstem 'LFPs/'], [], this_output_folder_tail,this_cond)
            fprintf('\n\nLFP DCM modelling complete for run number %d,\n\n',find(strcmp(subjcondpair{this_one,1},all_names)));
        catch
            try %Try again
                Preprocessing_mainfunction('extDCM',this_input_fname{thismeg},p,[pathstem 'LFPs/'], [], this_output_folder_tail,this_cond)
            catch
                fprintf('\n\nLFP DCM modelling failed for run number %d,\n\n',find(strcmp(subjcondpair{this_one,1},all_names)));
            end
        end
    end
end
    
subjcondpair = secondwavedata;
p.subjcntforcondition = 1;
p.multilevel = 1;
parfor this_one = 1:size(subjcondpair,1)
    this_input_fname = {['b8LFP_s_' time_wind_path{wind_cnt} '_' inv_meth{p.inv_cnt} '_' prefix subjcondpair{this_one,1} '.mat']};
    this_output_folder_tail = [Participant{find(strcmp(subjcondpair{this_one,1},all_names))}.diag '/']
    this_cond = find(strcmp(subjcondpair{this_one,2}, p.conditions));
    pause(mod(this_one*30,90)); %Introduce a pause to stagger the workers - otherwise sometimes the pool fails if trying to read or write simultaneously 
    for thismeg = 1:length(this_input_fname)
        try
            Preprocessing_mainfunction('extDCM',this_input_fname{thismeg},p,[pathstem 'LFPs/'], [], this_output_folder_tail,this_cond)
            fprintf('\n\nLFP DCM modelling complete for run number %d,\n\n',find(strcmp(subjcondpair{this_one,1},all_names)));
        catch
            try %Try again
                Preprocessing_mainfunction('extDCM',this_input_fname{thismeg},p,[pathstem 'LFPs/'], [], this_output_folder_tail,this_cond)
            catch
                fprintf('\n\nLFP DCM modelling failed for run number %d,\n\n',find(strcmp(subjcondpair{this_one,1},all_names)));
            end
        end
    end
end

delete(gcp)

%% Now do a first level PEB on the extDCM data
dirname_DCM = '/imaging/tc02/Holly_MMN/extDCMs_hE6/';
filestem = 'b8LFP_s_-100_500_LOR_fmbraedfffM';
conditions = {'STD','DVT','location','intensity','duration','gap','frequency'};
all_combinations = combvec(unique(p.group)',1:length(conditions));
if length(all_combinations) < 48
    Poolinfo = cbupool(length(all_combinations),'--mem-per-cpu=16G --time=167:00:00 --exclude=node-i[01-15]');
else
    Poolinfo = cbupool(48,'--mem-per-cpu=16G --time=167:00:00 --exclude=node-i[01-15]');
end
parpool(Poolinfo,Poolinfo.NumWorkers,'SpmdEnabled',false);
parfor this_comb = 1:length(all_combinations)
%for this_comb = 1:length(all_combinations) %falls over in parallel due to tmp.mat and unpredictable cd behaviour - needs fixing for bigger datasets
    k = all_combinations(1,this_comb)
    c = all_combinations(2,this_comb)
    extDCM_firstlevel_PEB(dirname_DCM,filestem,conditions(c),k,p,all_names)
end
rmdir([dirname_DCM 'PEB_firstlevel' filesep 'tempdir_*'])
delete(gcp)

%% Now do a second level PEB on the extDCM data
Poolinfo = cbupool(2*length(unique(p.group)'),'--mem-per-cpu=48G --time=167:00:00 --exclude=node-i[01-15]');
parpool(Poolinfo,Poolinfo.NumWorkers,'SpmdEnabled',false);
conditions = {'STD','DVT'}; %Tolerance failure if all conditions included
extDCM_secondlevel_PEB(dirname_DCM,filestem,conditions,unique(p.group)',p,all_names)
delete(gcp)

%% Now plot the whole scalp ERPs for sanity check
for todonumber = 1:nsubj
        try
        Participant{todonumber}.name = Participant{todonumber}.namepostmerge;
        end
end
prefix = 'PfmbraedfffM';
%plot_all_LFPs(Participant,pathstem,p,prefix)
% for i = biomarker_positive_mci_indices
%     Participant{i}.diag = 'MCI_pos';
% end
% for i = biomarker_negative_mci_indices
%     Participant{i}.diag = 'MCI_neg';
% end
% for i = biomarker_unknown_mci_indices
%     Participant{i}.diag = 'MCI_unk';
% end
% for i = biomarker_unknown_AD_indices
%     Participant{i}.diag = 'AD_unk';
% end
thesediagnoses = {'Control','ADMCI','nfvppa','pca','bvFTD'};
plot_ERP_bytype(Participant,pathstem,p,prefix,thesediagnoses)
quantify_MMN_ERP(Participant,pathstem,p,prefix,thesediagnoses)
%% Work in progress - quantitative analysis on the LFP data after Hughes et al. 2013

for todonumber = 1:nsubj
        try
        Participant{todonumber}.name = Participant{todonumber}.namepostmerge;
        end
end
prefix = 'fmbraedfffM';
val = 2; %for LORETA
%val = 1 %for IID
p.time_wind_path = time_wind_path;
p.wind_cnt = wind_cnt;
p.inv_meth = inv_meth;
p.inv_cnt = val;
baselined = 1;
plot_all_LFPs(Participant,pathstem,p,prefix)
quantify_MMN_LFP(Participant,pathstem,p,prefix,baselined)



% 
% prefix = 'PfmbraedfffM';
% secondlevelcomplete = zeros(1,1);
% this_output_folder_tail = {};
% p.mod = {'LFP'};
% for todonumber = 1:nsubj
%     if iscell(Participant{todonumber}.name)
%         this_output_folder_tail{todonumber}  = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.namepostmerge '/'];
%     else
%         this_output_folder_tail{todonumber}  = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/'];
%     end
% end
% p.all_conditions = p.conditions;
% p.conditions = {'STD','DVT'};
% try
%     tc_batch_SPM_LFP(prefix,this_output_folder_tail,pathstem,p);
%     secondlevelcomplete(1) = 1;
% catch
%     secondlevelcomplete(1) = 0;
% end
% p.conditions = {'STD','location','intensity','duration','gap','frequency'};
% try
%     tc_batch_SPM_LFP(prefix,this_output_folder_tail,pathstem,p);
%     secondlevelcomplete(1) = 1;
% catch
%     secondlevelcomplete(1) = 0;
% end
% p.conditions = p.all_conditions;

%% Bespoke connectivity test scripts
% TomHollyTA1
% compare_coherence_connectivity
