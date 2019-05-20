function maxfilter_this_participant(rawdatapath,subjfolder,denoisedfilename,startagain)
% e.g. maxfilter_this_participant(thesepaths,subjfolder,denoisedfilename)

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

clear all
clc
cd /imaging/hp02/pnfa_mmn/preprocessed

%% Add paths

addpath(genpath('/imaging/hp02/spm12b'));
addpath(genpath('/imaging/local/software/mne'));
addpath('/imaging/hp02/finger_tapping08/analysis_spm/new_functions');
addpath('/imaging/hp02/pnfa_mmn/');

%% Enter Subjects - Change Accordingly:
% script that holds all the participant's maxfilter locations and folder
% names
participant_folder_structure;
MF_path = '/imaging/hp02/pnfa_mmn/maxfilter/';

for ss = 1:length(Participant)
    mf_files{ss} = [MF_path, Participant{ss}.folder, '/', Participant{ss}.MF];
end


%% open matlabpool if required
ParType = 0;  % Fun on Login machines (not generally advised!)
%ParType = 1;   % Run maxfilter call on Compute machines using spmd (faster)

if ParType
cbupool(90);
end
%%
for ss =1: length(mf_files) 
    
    cd /imaging/hp02/pnfa_mmn/preprocessed/For_Thomas_dvts_sep/

    
    subjfolder = Participant{ss}.name;
        
    if ~exist( [Participant{ss}.folder, '/', subjfolder], 'dir'); mkdir([Participant{ss}.folder, '/', subjfolder]); end
    cd([Participant{ss}.folder, '/', subjfolder])
    
    fprintf(1, 'Subject: %s\n', subjfolder);
    
    dat_wd = '/imaging/hp02/pnfa_mmn/maxfilter/';
    
    rawfile  = fullfile(mf_files{ss});
    
    
    fprintf(1, 'Processing %s\n', rawfile);

    
    %%        Convert:
    S = [];
    S.dataset = rawfile;
    S.channels = 'all';
    S.outfile  = sprintf('%s', subjfolder); %_%s',filestem,maxflag)
    
    if ~exist(sprintf('f%s.mat',subjfolder), 'file')
        D = spm_eeg_convert(S);
    else
        D = spm_eeg_load(sprintf('%s.mat',subjfolder));
    end
    
    blk1_fid = D.fiducials;