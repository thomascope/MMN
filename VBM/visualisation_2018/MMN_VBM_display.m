%%
%% Display source reconstructions and bar charts for each group and condition

cfg.plots = [1:2];
cfg.symmetricity = 'symmetrical';
% cfg.normalise = 1;
% cfg.threshold = [5 40];
cfg.inflate = 10;

addpath([pwd '/ojwoodford-export_fig-216b30e'])

cfg.normalise = 0;

cfg.threshold = [4.7465 11]; %p=0.05 FWE

jp_spm8_surfacerender2_version_tc(['/imaging/tc02/Holly_MMN/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0001.nii'],'jet',cfg)
savepath = ['./VBM/Con-PCA_05FWE'];
eval(['export_fig ' savepath '.png -transparent -m2.5'])
close all

jp_spm8_surfacerender2_version_tc(['/imaging/tc02/Holly_MMN/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0002.nii'],'jet',cfg)
savepath = ['./VBM/Con-bvFTD_05FWE'];
eval(['export_fig ' savepath '.png -transparent -m2.5'])
close all

jp_spm8_surfacerender2_version_tc(['/imaging/tc02/Holly_MMN/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0003.nii'],'jet',cfg)
savepath = ['./VBM/Con-nfvPPA_05FWE'];
eval(['export_fig ' savepath '.png -transparent -m2.5'])
close all

jp_spm8_surfacerender2_version_tc(['/imaging/tc02/Holly_MMN/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0007.nii'],'jet',cfg)
savepath = ['./VBM/PCA-bvFTD_05FWE'];
eval(['export_fig ' savepath '.png -transparent -m2.5'])
close all

jp_spm8_surfacerender2_version_tc(['/imaging/tc02/Holly_MMN/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0009.nii'],'jet',cfg)
savepath = ['./VBM/PCA-nfvPPA_05FWE'];
eval(['export_fig ' savepath '.png -transparent -m2.5'])
close all

jp_spm8_surfacerender2_version_tc(['/imaging/tc02/Holly_MMN/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0011.nii'],'jet',cfg)
savepath = ['./VBM/bvFTD-nfvPPA_05FWE'];
eval(['export_fig ' savepath '.png -transparent -m2.5'])
close all

cfg.threshold = [3.19 11]; %p=0.001 

jp_spm8_surfacerender2_version_tc(['/imaging/tc02/Holly_MMN/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0001.nii'],'jet',cfg)
savepath = ['./VBM/Con-PCA_001unc'];
eval(['export_fig ' savepath '.png -transparent -m2.5'])
close all

jp_spm8_surfacerender2_version_tc(['/imaging/tc02/Holly_MMN/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0002.nii'],'jet',cfg)
savepath = ['./VBM/Con-bvFTD_001unc'];
eval(['export_fig ' savepath '.png -transparent -m2.5'])
close all

jp_spm8_surfacerender2_version_tc(['/imaging/tc02/Holly_MMN/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0003.nii'],'jet',cfg)
savepath = ['./VBM/Con-nfvPPA_001unc'];
eval(['export_fig ' savepath '.png -transparent -m2.5'])
close all

jp_spm8_surfacerender2_version_tc(['/imaging/tc02/Holly_MMN/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0007.nii'],'jet',cfg)
savepath = ['./VBM/PCA-bvFTD_001unc'];
eval(['export_fig ' savepath '.png -transparent -m2.5'])
close all

jp_spm8_surfacerender2_version_tc(['/imaging/tc02/Holly_MMN/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0009.nii'],'jet',cfg)
savepath = ['./VBM/PCA-nfvPPA_001unc'];
eval(['export_fig ' savepath '.png -transparent -m2.5'])
close all

jp_spm8_surfacerender2_version_tc(['/imaging/tc02/Holly_MMN/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0011.nii'],'jet',cfg)
savepath = ['./VBM/bvFTD-nfvPPA_001unc'];
eval(['export_fig ' savepath '.png -transparent -m2.5'])
close all