function preproc_this_participant(subjfolder,rawfilenames,startagain,p)
% e.g. preproc_this_participant(ss,[pathstem Participant{ss}.groupfolder '/' Participant{ss}.name '/M' Participant{ss}.name '.mat'])

%% Add paths

% addpath(genpath('/imaging/rowe/archive/users/hp02/spm12b'));
addpath(genpath('/imaging/local/software/mne'));
addpath('/imaging/rowe/archive/users/hp02/finger_tapping08/analysis_spm/new_functions');
addpath('/imaging/rowe/archive/users/hp02/pnfa_mmn/');

%% Enter Subjects - Change Accordingly:
% script that holds all the participant's maxfilter locations and folder
% names
cd(subjfolder)

% if startagain
%     delete f*
%     delete d*
%     delete e*
%     delete P*
%     delete r*
%     delete w*
%     delete s*
%     delete mr*
%     delete a*
% end
% if startagain
%     delete Pwfmra*
%     delete fmra*
%     delete mra*
%     delete wfmra*
%     delete s*
%     delete SPM*
% end
if startagain
    delete Pfmbra*
    delete Pwfmbra*
    delete fmbra*
    delete mbra*
    delete wfmbra*
    delete s*
    delete SPM*
    delete tf*
    delete bra*
end

if ~iscell(rawfilenames)
    rawfilenames = {rawfilenames};
end

for i = 1:length(rawfilenames)
    clear D
    D.fname = ['M' rawfilenames{i}];
    
    %% Notch filter
    S = [];
    S.D = D.fname;
    S.type = 'butterworth';
    S.band = 'stop';
    S.freq = [49 51];
    S.dir = 'twopass';
    S.order = 5;
    S.prefix = 'f';
    
    if ~exist(sprintf('ffM%s.mat',rawfilenames{i}), 'file')
        D = spm_eeg_filter(S);
    else
        D = spm_eeg_load(sprintf('fM%s.mat',rawfilenames{i}));
    end
    
    
    
    %% Low pass filter
    S = [];
    S.D = D.fname;
    S.type = 'butterworth';
    S.band = 'low';
    S.freq = 100;
    S.dir = 'twopass';
    S.order = 5;
    S.prefix = 'f';
    
    if ~exist(sprintf('fffM%s.mat',rawfilenames{i}), 'file')
        D = spm_eeg_filter(S);
    else
        D = spm_eeg_load(sprintf('ffM%s.mat',rawfilenames{i}));
    end%
    
    
    %% high pass filter
    S = [];
    S.D = D.fname;
    S.type = 'butterworth';
    S.band = 'high';
    S.freq = 0.5;
    S.dir = 'twopass';
    S.order = 5;
    S.prefix = 'f';
    
    if ~exist(sprintf('dfffM%s.mat',rawfilenames{i}), 'file')
        D = spm_eeg_filter(S);
    else
        D = spm_eeg_load(sprintf('fffM%s.mat',rawfilenames{i}));
    end
    %% Holly added: Downsampling:
    
    S = [];
    S.D = D.fname;
    S.fsample_new = 250;
    %D = spm_eeg_downsample(S);
    
    if ~exist(sprintf('edfffM%s.mat',rawfilenames{i}), 'file')
        D = spm_eeg_downsample(S);
    else
        D = spm_eeg_load(sprintf('dfffM%s.mat',rawfilenames{i}));
    end
    
    
    %% Epoching:
    disp('epoching')
    %spm('defaults', 'eeg');
    
    S = [];
    S.D = D.fname;
    S.timewin = p.timewin;
    S.trialdef(1).conditionlabel = 'STD';
    S.trialdef(1).eventvalue = 11; % select the MEG/EEG trigger codes which correspond to the stimuli/response you want to epoch around
    
    S.trialdef(2).conditionlabel = 'DVT';
    S.trialdef(2).eventvalue = [1,2,3,4,5,6,7,8];
    
    S.trialdef(3).conditionlabel = 'location';
    S.trialdef(3).eventvalue = [7,8];
    
    S.trialdef(4).conditionlabel = 'intensity';
    S.trialdef(4).eventvalue = [5,6];
    
    S.trialdef(5).conditionlabel = 'duration';
    S.trialdef(5).eventvalue = 1;
    
    S.trialdef(6).conditionlabel = 'gap';
    S.trialdef(6).eventvalue = 4;
    
    S.trialdef(7).conditionlabel = 'frequency';
    S.trialdef(7).eventvalue = [2,3];
    
    S.trialdef(8).conditionlabel = 'location_L';
    S.trialdef(8).eventvalue = 7;
    
    S.trialdef(9).conditionlabel = 'frequency_high';
    S.trialdef(9).eventvalue = 3;
    
    S.trialdef(10).conditionlabel = 'intensity_high';
    S.trialdef(10).eventvalue = 6;
    
    S.trialdef(11).conditionlabel = 'location_R';
    S.trialdef(11).eventvalue = 8;
    
    S.trialdef(12).conditionlabel = 'intensity_low';
    S.trialdef(12).eventvalue = 5;
    
    S.trialdef(13).conditionlabel = 'frequency_low';
    S.trialdef(13).eventvalue = 2;
    
    for c=1:length(S.trialdef)
        S.trialdef(c).eventtype = 'STI101_up'; % look for rising signal in trigger channel for event times
        S.trialdef(c).trlshift = 0;
    end
    
    S.bc = 1;
    S.prefix = 'e';
    S.eventpadding = 0;
    %D = spm_eeg_epochs(S);
    
    if ~exist(sprintf('aedfffM%s.mat',rawfilenames{i}), 'file')
        D = spm_eeg_epochs(S);
    else
        D = spm_eeg_load(sprintf('edfffM%s.mat',rawfilenames{i}));
    end
    
    %% Artifact rejection
    
    S = [];
    S.D = D.fname;
    S.mode = 'reject';
    S.badchanthresh = 0.2;
    
    %
    S.methods(1).channels = {'EOG'};
    S.methods(1).fun = 'peak2peak';
    %S.methods(1).fun = 'threshchan';
    S.methods(1).settings.threshold = 200;
    
    
    S.methods(end).channels = {'MEGPLANAR'};
    S.methods(end+1).fun = 'peak2peak';
    S.methods(end).settings.threshold = 200; %900ft%meg_thr(ss);
    %
    %         S.methods(end).channels = {'EEG'};
    %         S.methods(end+1).fun = 'peak2peak';
    %         S.methods(end).settings.threshold = 200;
    
    %D = spm_eeg_artefact(S);
    
    if ~exist(sprintf('raedfffM%s.mat',rawfilenames{i}), 'file')
        D = spm_eeg_artefact(S);
    else
        D = spm_eeg_load(sprintf('aedfffM%s.mat',rawfilenames{i}))
    end
    
    
    
    %nbadchans(ss) = length(D.badchannels);
    %nrej = 0;
    
    disp('done');
    
    %% Remove Bad channels observed in experiment and maxfilter
    
    
    %D = badchannels(D, indchannel(D,badchannels{ss,1}), 1);
    %parsave(D.fname, D);
    
    %% Remove bad trials:
    S = [];
    S.D = D.fname;
    S.prefix = 'r';
    
    %D = spm_eeg_remove_bad_trials(S);
    if ~exist(sprintf('braedfffM%s.mat',rawfilenames{i}), 'file')
        D = spm_eeg_remove_bad_trials(S);
    else
        D = spm_eeg_load(sprintf('raedfffM%s.mat',rawfilenames{i}));
    end
    
    %% Baseline correct
    S = [];
    S.timewin = [p.preBase p.postBase];
    S.D = D.fname;
    S.prefix = 'b';
    D = spm_eeg_bc(S);    
end

%% Merge
if length(rawfilenames)>1 && ~exist(sprintf('mbraedfffM%s.mat',rawfilenames{1}(1:end-2)), 'file')
   basename = D.fname;
   basename = basename(1:end-6); % Assume fewer than 10 files to merge and naming convention _x
   S = [];
   S.D = {};
   for i = 1:length(rawfilenames)
       S.D{i} = [basename '_' num2str(i) '.mat'];
   end
   S.D = char(S.D);
   D = spm_eeg_merge(S);
   postmerge_fname = D.fname;
   S.D = postmerge_fname;
   S.outfile = [basename '.mat'];
   D = spm_eeg_copy(S);
   delete([postmerge_fname(1:end-4) '*'])
end

if length(rawfilenames)>1
    rawfilenames = {rawfilenames{1}(1:end-2)};
end

%% Robust averaging
S = [];
S.D = D.fname;
S.robust.ks = 3;
S.robust.bycondition = true;
S.robust.savew = false;
S.robust.removebad = false;
S.circularise = false;
S.prefix = 'm';
% D = spm_eeg_average(S);

if ~exist(sprintf('fmbraedfffM%s.mat',rawfilenames{1}), 'file')
    D = spm_eeg_average(S);
else
    D = spm_eeg_load(sprintf('mbraedfffM%s.mat',rawfilenames{1}));
end


%% Low pass filter (again for robust averaging)
S = [];
S.D = D.fname;
S.type = 'butterworth';
S.band = 'low';
S.freq = 40;
S.dir = 'twopass';
S.order = 5;
S.prefix = 'f';
%         D = spm_eeg_filter(S);

if ~exist(sprintf('wfmbraedfffM%s.mat',rawfilenames{1}), 'file')
    D = spm_eeg_filter(S);
else
    D = spm_eeg_load(sprintf('fmbraedfffM%s.mat',rawfilenames{1}));
end

%% Compute Contrasts
S = [];
S.D = D.fname;
%S.c = [1 -1];
S.c = p.contrast_weights;
S.label = p.contrast_labels;

%S.label = {'std-dvt'};



S.weighted = 1;
S.prefix = 'w';
%     D = spm_eeg_contrast(S);
if ~exist(sprintf('PwfmbraedfffM%s.mat',rawfilenames{1}), 'file')
    D = spm_eeg_contrast(S);
else
    D = spm_eeg_load(sprintf('wfmbraedfffM%s.mat',rawfilenames{1}));
end

%% Combine Planars for MMN

S = [];
S.D = D.fname;
S.mode = 'replace';
D = spm_eeg_combineplanar(S);


%% Combine Planars for std and dvt

S = [];
S.D = sprintf('fmbraedfffM%s.mat',rawfilenames{1});
S.mode = 'replace';
D = spm_eeg_combineplanar(S);

%% Delete files to keep space usage down

% delete fM*
% delete ff*
% delete df*
% delete ed*
% delete ae*
% delete mr*
% delete df*


cd ../

return