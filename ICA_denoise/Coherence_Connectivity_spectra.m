function [frequency_spectra, frequencies, regions] = Coherence_Connectivity_spectra(Participant,pathstem,p,prefix)
%method = 'granger';
%method = 'coh';

rmpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/'))
%addpath /imaging/local/software/spm_cbu_svn/releases/spm12_fil_r6906
addpath /group/language/data/thomascope/spm12_fil_r6906/
addpath /group/language/data/thomascope/MMN/Coherence_Connectivity_Analysis
%spm eeg

thisdir = pwd;

addpath('/group/language/data/thomascope/vespa/SPM12version/Standalone preprocessing pipeline/')

pathstem_LFPs = '/imaging/mlr/users/tc02/Holly_MMN/ICA_denoise/LFPs';

outdir = ['/imaging/mlr/users/tc02/Holly_MMN/Coherence_Connectivity_Integrated_' p.inv_meth{p.inv_cnt} '/'];

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
    megpath{ss} = [pathstem Participant{ss}.groupfolder '/' Participant{ss}.namepostmerge '/' 's_' p.time_wind_path{p.wind_cnt} '_' p.inv_meth{p.inv_cnt} '_' prefix Participant{ss}.namepostmerge '.mat'];
    diagnosis{ss} = Participant{ss}.diag;
    
    [f1,f2,f3] = fileparts(megpath{ss});
    if strcmp(diagnosis{ss},'AD_MCI_pos')
        fn{ss} = sprintf('%s/%s/%dLFP_%s%s',[pathstem 'LFPs'],'MCI', length(Sname), f2, f3);
    else
        fn{ss} = sprintf('%s/%s/%dLFP_%s%s',[pathstem 'LFPs'],diagnosis{ss}, length(Sname), f2, f3);
    end
    if ~exist(fn{ss},'file')
        if exist([fn{ss}(1:end-4) '_1.mat'])
            movefile([fn{ss}(1:end-4) '_1.mat'],fn{ss})
        else
            error(['File ' fn{ss} ' not found'])
        end
    end
    
end

fft_method = 'mtmfft'; % 'wavelet' for morlet; can leave blank for multitaper.

[groupstodo,~, all_subjs] = unique(diagnosis,'stable');

decompositionworkedcorrectly = zeros(1,length(Participant));
for subj = 1:length(Participant)
    warning('off','all')
    
    [frequency_spectra(subj,:,:,:,:), frequencies, regions] = MMN_coherence_spectra(fn{subj},[],[],subj,start_times,end_times,fft_method,method);
end