function preproc_this_participant(subjfolder,denoisedfilename)
% e.g. preproc_this_participant(ss,[pathstem Participant{ss}.groupfolder '/' Participant{ss}.name '/M' Participant{ss}.name '.mat'])

%% Add paths

addpath(genpath('/imaging/hp02/spm12b'));
addpath(genpath('/imaging/local/software/mne'));
addpath('/imaging/hp02/finger_tapping08/analysis_spm/new_functions');
addpath('/imaging/hp02/pnfa_mmn/');

%% Enter Subjects - Change Accordingly:
% script that holds all the participant's maxfilter locations and folder
% names
cd(subjfolder)

D.fname = denoisedfilename;

%% Notch filter
S = [];
S.D = D.fname;
S.type = 'butterworth';
S.band = 'stop';
S.freq = [49 51];
S.dir = 'twopass';
S.order = 5;
S.prefix = 'f';

if ~exist(sprintf('ff%s.mat',denoisedfilename), 'file')
    D = spm_eeg_filter(S);
else
    D = spm_eeg_load(sprintf('f%s.mat',denoisedfilename));
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

if ~exist(sprintf('fff%s.mat',denoisedfilename), 'file')
    D = spm_eeg_filter(S);
else
    D = spm_eeg_load(sprintf('ff%s.mat',denoisedfilename));
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

if ~exist(sprintf('dfff%s.mat',denoisedfilename), 'file')
    D = spm_eeg_filter(S);
else
    D = spm_eeg_load(sprintf('fff%s.mat',denoisedfilename));
end
%% Holly added: Downsampling:

S = [];
S.D = D.fname;
S.fsample_new = 250;
%D = spm_eeg_downsample(S);

if ~exist(sprintf('edfff%s.mat',denoisedfilename), 'file')
    D = spm_eeg_downsample(S);
else
    D = spm_eeg_load(sprintf('dfff%s.mat',denoisedfilename));
end


%% Epoching:
disp('epoching')
%spm('defaults', 'eeg');

S = [];
S.D = D.fname;
S.timewin = [-100 500];
S.trialdef(1).conditionlabel = 'STD';
S.trialdef(1).eventvalue = 11; % select the MEG/EEG trigger codes which correspond to the stimuli/response you want to epoch around

S.trialdef(2).conditionlabel = 'DVT';
S.trialdef(2).eventvalue = [1,2,3,4,5,6,7,8];


for c=1:2
    S.trialdef(c).eventtype = 'STI101_up'; % look for rising signal in trigger channel for event times
    S.trialdef(c).trlshift = 0;
end

S.bc = 1;
S.prefix = 'e';
S.eventpadding = 0;
%D = spm_eeg_epochs(S);

if ~exist(sprintf('aedfff%s.mat',denoisedfilename), 'file')
    D = spm_eeg_epochs(S);
else
    D = spm_eeg_load(sprintf('edfff%s.mat',denoisedfilename));
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

if ~exist(sprintf('raedfff%s.mat',denoisedfilename), 'file')
    D = spm_eeg_artefact(S);
else
    D = spm_eeg_load(sprintf('aedfff%s.mat',denoisedfilename))
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
if ~exist(sprintf('mraedfff%s.mat',denoisedfilename), 'file')
    D = spm_eeg_remove_bad_trials(S);
else
    D = spm_eeg_load(sprintf('raedfff%s.mat',denoisedfilename));
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

if ~exist(sprintf('fmraedfff%s.mat',denoisedfilename), 'file')
    D = spm_eeg_average(S);
else
    D = spm_eeg_load(sprintf('mraedfff%s.mat',denoisedfilename));
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

if ~exist(sprintf('wfmraedfff%s.mat',denoisedfilename), 'file')
    D = spm_eeg_filter(S);
else
    D = spm_eeg_load(sprintf('fmraedfff%s.mat',denoisedfilename));
end

%% Compute Contrasts
S = [];
S.D = D.fname;
S.c = [1 -1];

S.label = {'std-dvt'};



S.weighted = 1;
S.prefix = 'w';
%     D = spm_eeg_contrast(S);
if ~exist(sprintf('wfmraedfff%s.mat',denoisedfilename), 'file')
    D = spm_eeg_contrast(S);
else
    D = spm_eeg_load(sprintf('wfmraedfff%s.mat',denoisedfilename));
end

%% Combine Planars for MMN

S = [];
S.D = D.fname;
S.mode = 'replace';
D = spm_eeg_combineplanar(S);


%% Combine Planars for std and dvt

S = [];
S.D = sprintf('fmraedfff%s.mat',denoisedfilename);
S.mode = 'replace';
D = spm_eeg_combineplanar(S);

%% Delete files to keep space usage down

delete Mm*
delete fM*
delete ff*
delete df*
delete ed*
delete ae*
delete mr*
delete df*


cd ../

return