%%
%% Display source reconstructions and bar charts for each group and condition

cfg.plots = [1:2];
cfg.symmetricity = 'symmetrical';
% cfg.normalise = 1;
% cfg.threshold = [5 40];
cfg.inflate = 10;

addpath([pwd '/ojwoodford-export_fig-216b30e'])

cfg.normalise = 0;

conditions = {'STD','DVT','Loc','Int','Dur','Gap','Freq'};

mkdir('./VBM_bycond')

cfg.threshold = [3.19 11]; %p=0.001 

this_file = ['/imaging/tc02/Holly_MMN/ICA_Denoise/VBM_STD_amplitude_left A1_M100/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0001.nii'];
jp_spm8_surfacerender2_version_tc(this_file,'jet',cfg)
savepath = ['./VBM_bycond/standard_M100_001unc'];
eval(['export_fig ' savepath '.png -transparent -m2.5'])
close all


for i = 2:length(conditions)
    this_file = ['/imaging/tc02/Holly_MMN/ICA_Denoise/VBM_MMN_amplitude_left A1_' conditions{i} '/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0001.nii'];
    jp_spm8_surfacerender2_version_tc(this_file,'jet',cfg)
    savepath = ['./VBM_bycond/' conditions{i} '_amplitude_001unc'];
    eval(['export_fig ' savepath '.png -transparent -m2.5'])
    close all
    
    this_file = ['/imaging/tc02/Holly_MMN/ICA_Denoise/VBM_MMN_latency_left A1_' conditions{i} '/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0001.nii'];
    jp_spm8_surfacerender2_version_tc(this_file,'jet',cfg)
    savepath = ['./VBM_bycond/' conditions{i} '_latency_001unc'];
    eval(['export_fig ' savepath '.png -transparent -m2.5'])
    close all

end

cfg.threshold = [4.7465 11]; %p=0.05 FWE

this_file = ['/imaging/tc02/Holly_MMN/ICA_Denoise/VBM_STD_amplitude_left A1_M100/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0001.nii'];
jp_spm8_surfacerender2_version_tc(this_file,'jet',cfg)
savepath = ['./VBM_bycond/standard_M100_05FWE'];
eval(['export_fig ' savepath '.png -transparent -m2.5'])
close all


for i = 2:length(conditions)
    this_file = ['/imaging/tc02/Holly_MMN/ICA_Denoise/VBM_MMN_amplitude_left A1_' conditions{i} '/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0001.nii'];
    jp_spm8_surfacerender2_version_tc(this_file,'jet',cfg)
    savepath = ['./VBM_bycond/' conditions{i} '_amplitude_05FWE'];
    eval(['export_fig ' savepath '.png -transparent -m2.5'])
    close all
    
    this_file = ['/imaging/tc02/Holly_MMN/ICA_Denoise/VBM_MMN_latency_left A1_' conditions{i} '/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0001.nii'];
    jp_spm8_surfacerender2_version_tc(this_file,'jet',cfg)
    savepath = ['./VBM_bycond/' conditions{i} '_latency_05FWE'];
    eval(['export_fig ' savepath '.png -transparent -m2.5'])
    close all

end

%% Now extract the t-scores from each ROI
Lpos = [-42, -22, 7;
-61, -32, 8;
-46, 20, 8;
-49, -38, 38;
46, -14, 8;
59, -25, 8;
46, 20, 8;
57, -38, 42];

this_file = ['/imaging/tc02/Holly_MMN/ICA_Denoise/VBM_STD_amplitude_left A1_M100/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0001.nii'];
mean_tscore = zeros(length(conditions),8);
for i = 1:size(Lpos,1)
this_tmatrix = spm_summarise(this_file,struct('def','sphere', 'spec',8, 'xyz',Lpos(i,:)'));
mean_tscore(1,i) = mean(nonzeros(this_tmatrix));
end


for j = 2:length(conditions)
    this_file = ['/imaging/tc02/Holly_MMN/ICA_Denoise/VBM_MMN_amplitude_left A1_' conditions{j} '/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0001.nii'];
    for i = 1:size(Lpos,1)
        this_tmatrix = spm_summarise(this_file,struct('def','sphere', 'spec',8, 'xyz',Lpos(i,:)'));
        mean_tscore(j,i) = mean(nonzeros(this_tmatrix));
    end

    
    %this_file = ['/imaging/tc02/Holly_MMN/ICA_Denoise/VBM_MMN_latency_left A1_' conditions{i} '/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0001.nii'];

end



% 
% jp_spm8_surfacerender2_version_tc(['/imaging/tc02/Holly_MMN/ICA_Denoise/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0001.nii'],'jet',cfg)
% savepath = ['./VBM/Con-PCA_001unc'];
% eval(['export_fig ' savepath '.png -transparent -m2.5'])
% close all
% 
% jp_spm8_surfacerender2_version_tc(['/imaging/tc02/Holly_MMN/ICA_Denoise/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0002.nii'],'jet',cfg)
% savepath = ['./VBM/Con-bvFTD_001unc'];
% eval(['export_fig ' savepath '.png -transparent -m2.5'])
% close all
% 
% jp_spm8_surfacerender2_version_tc(['/imaging/tc02/Holly_MMN/ICA_Denoise/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0003.nii'],'jet',cfg)
% savepath = ['./VBM/Con-nfvPPA_001unc'];
% eval(['export_fig ' savepath '.png -transparent -m2.5'])
% close all
% 
% jp_spm8_surfacerender2_version_tc(['/imaging/tc02/Holly_MMN/ICA_Denoise/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0007.nii'],'jet',cfg)
% savepath = ['./VBM/PCA-bvFTD_001unc'];
% eval(['export_fig ' savepath '.png -transparent -m2.5'])
% close all
% 
% jp_spm8_surfacerender2_version_tc(['/imaging/tc02/Holly_MMN/ICA_Denoise/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0009.nii'],'jet',cfg)
% savepath = ['./VBM/PCA-nfvPPA_001unc'];
% eval(['export_fig ' savepath '.png -transparent -m2.5'])
% close all
% 
% jp_spm8_surfacerender2_version_tc(['/imaging/tc02/Holly_MMN/ICA_Denoise/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0011.nii'],'jet',cfg)
% savepath = ['./VBM/bvFTD-nfvPPA_001unc'];
% eval(['export_fig ' savepath '.png -transparent -m2.5'])
% close all
% 
% jp_spm8_surfacerender2_version_tc(['/imaging/tc02/vespa/scans/PNFA_VBM/tom/Stats/ERF_factorial_full_group_vbm_TIVnormalised_agecovaried_smoothedmask/spmT_0001.nii'],'jet',cfg)
% savepath = ['./VBM/Con-nfvPPA_VESPA_001unc'];
% eval(['export_fig ' savepath '.png -transparent -m2.5'])
% close all
% 
% jp_spm8_surfacerender2_version_tc(['/imaging/tc02/vespa/scans/PNFA_VBM/tom/Stats/ERF_factorial_full_group_vbm_TIVnormalised_agecovaried_smoothedmask/spmT_0002.nii'],'jet',cfg)
% savepath = ['./VBM/Con-nfvPPA_VESPA_001unc_rev'];
% eval(['export_fig ' savepath '.png -transparent -m2.5'])
% close all
% 
% jp_spm8_surfacerender2_version_tc(['/imaging/tc02/Holly_MMN/ICA_Denoise/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0013.nii'],'jet',cfg)
% savepath = ['./VBM/Con-ADMCI_001unc'];
% eval(['export_fig ' savepath '.png -transparent -m2.5'])
% close all
% 
% jp_spm8_surfacerender2_version_tc(['/imaging/tc02/Holly_MMN/ICA_Denoise/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0015.nii'],'jet',cfg)
% savepath = ['./VBM/PCA-ADMCI_001unc'];
% eval(['export_fig ' savepath '.png -transparent -m2.5'])
% close all
% 
% jp_spm8_surfacerender2_version_tc(['/imaging/tc02/Holly_MMN/ICA_Denoise/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0017.nii'],'jet',cfg)
% savepath = ['./VBM/bvFTD-ADMCI_001unc'];
% eval(['export_fig ' savepath '.png -transparent -m2.5'])
% close all
% 
% jp_spm8_surfacerender2_version_tc(['/imaging/tc02/Holly_MMN/ICA_Denoise/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0019.nii'],'jet',cfg)
% savepath = ['./VBM/nfvPPA-ADMCI_001unc'];
% eval(['export_fig ' savepath '.png -transparent -m2.5'])
% close all
