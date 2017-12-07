function Master_Script_secondfilter

%A script to calculate Imaginary Coherence and Granger Causailty from LFP
%data for the MMN project

makefiles = 0;
extractlfps = 0;

method = 'granger';
%method = 'coh';

start_times = 0;
end_times = 500;

thisdir = pwd;

addpath('/group/language/data/thomascope/vespa/SPM12version/Standalone preprocessing pipeline/')

%pathstem = '/imaging/hp02/pnfa_mmn/forward_modelling/LFPs_TA';
pathstem_evoked = '/imaging/tc02/Holly_MMN/LFPs/';
pathstem_inv = '/imaging/hp02/pnfa_mmn/forward_modelling/MSP/';
pathstem_raw = '/imaging/hp02/pnfa_mmn/preprocessed/For_Thomas_dvts_sep/';
outdir = '/imaging/tc02/Holly_MMN/Coherence_Connectivity_secondfilter/';

fft_method = 'mtmfft'; % 'wavelet' for morlet; can leave blank for multitaper.

groupstodo = {'matched_HCs' 'pca' 'bvFTD' 'pnfa'};
dirnames_inv = {'matched_HCs' 'pca' 'ftd' 'vespa'};
all_subjs = [];

cd(pathstem_evoked)
filenames = {};

for i = 1:length(groupstodo)
    cd(groupstodo{i});
    thesefiles = dir('8LFP*mat');
    all_subjs = [all_subjs; [repmat(i,size(thesefiles)), (1:size(thesefiles,1))']];
    filenames = [filenames, thesefiles.name];
    cd(pathstem_evoked)
end

inversion_prefix = 'fmraedfff';

raw_prefix = 'raedffff';

val = 1; %1st inversion

if makefiles == 1
for s=1:size(all_subjs,1)
    %Find out where all the data are
    raw_suffix = strsplit(filenames{s},inversion_prefix);
    raw_suffix = raw_suffix{2}; %get the participant i.d.
    if exist([pathstem_inv dirnames_inv{all_subjs(s,1)} '/' raw_suffix(1:end-4) '/' inversion_prefix raw_suffix],'file');
        S.D = [pathstem_inv dirnames_inv{all_subjs(s,1)} '/' raw_suffix(1:end-4) '/' inversion_prefix raw_suffix];    
    else
        error([pathstem_inv dirnames_inv{all_subjs(s,1)} '/' raw_suffix(1:end-4)  ' does not contain the inversion file'])
    end
    if exist([pathstem_raw dirnames_inv{all_subjs(s,1)} '/' raw_suffix(1:end-4) '/' raw_prefix raw_suffix],'file')
        unaveragedfile = [pathstem_raw dirnames_inv{all_subjs(s,1)} '/' raw_suffix(1:end-4) '/' raw_prefix raw_suffix];
    else
        error([pathstem_inv dirnames_inv{all_subjs(s,1)} '/' raw_suffix(1:end-4)  ' does not contain the raw data file'])
    end
    ouputfile = [outdir groupstodo{all_subjs(s,1)} '/' raw_suffix(1:end-4) '/' raw_prefix raw_suffix];
    
    if ~exist([outdir groupstodo{all_subjs(s,1)} '/' raw_suffix(1:end-4) '/'],'dir')
        mkdir([outdir groupstodo{all_subjs(s,1)} '/' raw_suffix(1:end-4) '/'])
    end
    
    %Extract the inversion
    D = spm_eeg_load(S.D);
    inversion = D.inv;
    
    clear S D
    
    %Copy the raw data
    S.D = unaveragedfile;
    S.outfile = ouputfile;
    
    spm_eeg_copy(S)
    
    %Place the inversion in raw data file
    D = spm_eeg_load(S.outfile);
    D.inv = inversion;
    D.save

end
end



cd(thisdir)
s = matlabpool('size');
if s == 0 
this_pool = cbupool(size(all_subjs,1));
this_pool.ResourceTemplate = '-l nodes=^N^,mem=256GB,walltime=100:00:00';
matlabpool(this_pool)
end

if extractlfps == 1
    parfor subj = 1:size(all_subjs,1)
        raw_suffix = strsplit(filenames{subj},inversion_prefix);
        raw_suffix = raw_suffix{2}; %get the participant i.d.
        this_rawfile = [outdir groupstodo{all_subjs(subj,1)} '/' raw_suffix(1:end-4) '/' raw_prefix raw_suffix]
        
        D = spm_eeg_load(this_rawfile)
    D.inv{val}.source.XYZ = [-42, -22, 7;
        -61, -32, 8;
        -46, 20, 8;
        -49, -38, 38;
        46, -14, 8;
        59, -25, 8;
        46, 20, 8;
        57, -38, 42];
    D.inv{val}.source.label = {'left A1';
        'left STG';
        'left IFG';
        'left IPC';
        'right A1';
        'right STG';
        'right IFG';
        'right IPC'};
    D.inv{val}.source.fname = [outdir 'timeseries_for_coherence_' groupstodo{all_subjs(subj,1)} '_' num2str(all_subjs(subj,2))];
    D.inv{val}.source.type = 'trials';
    D.val = val;
    thisdir = pwd;
    cd([outdir groupstodo{all_subjs(subj,1)} '/' raw_suffix(1:end-4) '/']);
    spm_eeg_inv_extract(D)
    cd(thisdir);
   
    end
end

parfor subj = 1:size(all_subjs,1)
    warning('off','all')
    if strcmp(method,'coh')
        this_outdir = [outdir groupstodo{all_subjs(subj,1)} '/'];
    else
        this_outdir = [outdir method '/' groupstodo{all_subjs(subj,1)} '/'];
    end
    
    if exist([this_outdir 's' num2str(all_subjs(subj,:)) '_grangerdata_highfreq_averagesubtracted_100_' num2str(start_times) '_' num2str(end_times) '_z' num2str(1) '.mat'],'file')
        disp(['Found some existing files for subject s' num2str(all_subjs(subj,:))])
        existing_files = dir([this_outdir 's' num2str(all_subjs(subj,:)) '_grangerdata_highfreq_averagesubtracted_100_' num2str(start_times) '_' num2str(end_times) '_z*'])
        this_file_num = [];
        for i = 1:length(existing_files)
        temp_files = strsplit(existing_files(i).name,'z');
        temp_files = strsplit(temp_files{2},'.mat');
        this_file_num = [this_file_num str2num(temp_files{1})];
        end
        uptohere = max(this_file_num)+1;
        disp(['Resuming at file ' num2str(uptohere)])
        parallel_MMN_coherence_granger_resume([outdir 'timeseries_for_coherence_' groupstodo{all_subjs(subj,1)} '_' num2str(all_subjs(subj,2))],[],this_outdir,all_subjs(subj,:),start_times,end_times,fft_method,method,uptohere)
    else
        disp(['job started on worker ' num2str(subj)])
        parallel_MMN_coherence_granger([outdir 'timeseries_for_coherence_' groupstodo{all_subjs(subj,1)} '_' num2str(all_subjs(subj,2))],[],this_outdir,all_subjs(subj,:),start_times,end_times,fft_method,method)
    end
end

matlabpool close