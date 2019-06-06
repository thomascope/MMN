function decompositionworkedcorrectly = Coherence_Connectivity(Participant,pathstem,p,prefix)
%method = 'granger';
%method = 'coh';

rmpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/'))
%addpath /imaging/local/software/spm_cbu_svn/releases/spm12_fil_r6906
addpath /group/language/data/thomascope/spm12_fil_r6906/
addpath /group/language/data/thomascope/MMN/Coherence_Connectivity_Analysis
%spm eeg

thisdir = pwd;

addpath('/group/language/data/thomascope/vespa/SPM12version/Standalone preprocessing pipeline/')

pathstem_LFPs = '/imaging/tc02/Holly_MMN/ICA_denoise/LFPs';
outdir = '/imaging/tc02/Holly_MMN/Coherence_Connectivity_ICA/';

%A script to calculate Imaginary Coherence and Granger Causailty from LFP
%data for the MMN project
start_times = p.start_times;
end_times = p.end_times;

method = p.decompmethod;

%
Sname = {'left A1';
    'left STG';
    'left IFG';
    'left IPC';
    'right A1';
    'right STG';
    'right IFG';
    'right IPC'};

for ss = 1:length(Participant)
    megpath{ss} = [pathstem Participant{ss}.groupfolder '/' Participant{ss}.name '/' 's_' p.time_wind_path{p.wind_cnt} '_' p.inv_meth{p.inv_cnt} '_' prefix Participant{ss}.name '.mat'];
    diagnosis{ss} = Participant{ss}.diag;
    
    [f1,f2,f3] = fileparts(megpath{ss});
    fn{ss} = sprintf('%s/%s/%dLFP_%s%s',[pathstem 'LFPs'],diagnosis{ss}, length(Sname), f2, f3);
    
end

fft_method = 'mtmfft'; % 'wavelet' for morlet; can leave blank for multitaper.

[groupstodo,~, all_subjs] = unique(diagnosis,'stable');

decompositionworkedcorrectly = zeros(1,length(Participant));
parfor subj = 1:length(Participant)
    warning('off','all')
    if strcmp(method,'coh')
        this_outdir = [outdir 'crosshem/' method '/' groupstodo{all_subjs(subj,1)} '/'];
    else
        this_outdir = [outdir 'crosshem/' method '/' groupstodo{all_subjs(subj,1)} '/'];
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
        try
        parallel_MMN_coherence_granger_resume(fn{subj},[],this_outdir,subj,start_times,end_times,fft_method,method,uptohere)
        decompositionworkedcorrectly(subj) = 1;
        end
    else
        disp(['job started on worker ' num2str(subj)])
        try
        parallel_MMN_coherence_granger(fn{subj},[],this_outdir,subj,start_times,end_times,fft_method,method)
        decompositionworkedcorrectly(subj) = 1;
        end
    end
end