%%
%% Display source reconstructions and bar charts for each group and condition

cfg.plots = [1:2];
cfg.symmetricity = 'symmetrical';
% cfg.normalise = 1;
% cfg.threshold = [5 40];
cfg.inflate = 10;

addpath(genpath('/group/language/data/thomascope/MMN/ICA_denoise/VBM/visualisation_2020'))

cfg.normalise = 0;

cfg.threshold = [10.9368 31.73]; %p=0.0001 unc

if ~exist('./SPM','dir')
    mkdir('./SPM')
end

jp_spm8_surfacerender2_version_tc(['/imaging/mlr/users/tc02/Holly_MMN/ICA_denoise_longwindow/source_stats/-100_500/spmF_0002.nii'],'jet',cfg)
savepath = ['./SPM/MMN_All_groups'];
eval(['export_fig ' savepath '.png -transparent -m2.5'])
close all
