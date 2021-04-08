%%
%% Display source reconstructions and bar charts for each group and condition

cfg.plots = [1:2];
cfg.symmetricity = 'symmetrical';
% cfg.normalise = 1;
% cfg.threshold = [5 40];
cfg.inflate = 10;

addpath([pwd '/ojwoodford-export_fig-216b30e'])

cfg.normalise = 0;

mkdir('./VBM_connectivity')

cfg.threshold = [3.19 11]; %p=0.001 

Granger_dirs = dir('/imaging/mlr/users/tc02/Holly_MMN/ICA_Denoise/VBM_G*');
icoh_plv_dirs = dir('/imaging/mlr/users/tc02/Holly_MMN/ICA_Denoise/VBM_partial*');

Lpos = [-42, -22, 7;
-61, -32, 8;
-46, 20, 8;
-49, -38, 38;
46, -14, 8;
59, -25, 8;
46, 20, 8;
57, -38, 42];

mean_tscore = zeros(length(Granger_dirs) + length(icoh_plv_dirs),8);
max_tscore = zeros(length(Granger_dirs) + length(icoh_plv_dirs),8);

for i = 1:length(Granger_dirs)
    this_file = ['/imaging/mlr/users/tc02/Holly_MMN/ICA_Denoise/' Granger_dirs(i).name '/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0001.nii'];
    jp_spm8_surfacerender2_version_tc(this_file,'jet',cfg)
    savepath = ['./VBM_connectivity/' Granger_dirs(i).name '_001unc'];
    eval(['export_fig ' strrep(savepath, ' ', '_') '.png -transparent -m2.5'])
    close all
    for j = 1:size(Lpos,1)
        this_tmatrix = spm_summarise(this_file,struct('def','sphere', 'spec',8, 'xyz',Lpos(j,:)'));
        mean_tscore(i,j) = mean(nonzeros(this_tmatrix));
        max_tscore(i,j) = max(abs(nonzeros(this_tmatrix)));
    end
end

for i = 1:length(icoh_plv_dirs)
    this_file = ['/imaging/mlr/users/tc02/Holly_MMN/ICA_Denoise/' icoh_plv_dirs(i).name '/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0001.nii'];
    jp_spm8_surfacerender2_version_tc(this_file,'jet',cfg)
    savepath = ['./VBM_connectivity/' icoh_plv_dirs(i).name '_001unc'];
    eval(['export_fig ' strrep(savepath, ' ', '_') '.png -transparent -m2.5'])
    close all
    for j = 1:size(Lpos,1)
        this_tmatrix = spm_summarise(this_file,struct('def','sphere', 'spec',8, 'xyz',Lpos(j,:)'));
        mean_tscore(length(Granger_dirs)+i,j) = mean(nonzeros(this_tmatrix));
        max_tscore(length(Granger_dirs)+i,j) = max(abs(nonzeros(this_tmatrix)));
    end
end


% %% Now extract the t-scores from each ROI
% 
% 
% this_file = ['/imaging/mlr/users/tc02/Holly_MMN/ICA_Denoise/VBM_STD_amplitude_left A1_M100/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0001.nii'];
% mean_tscore = zeros(length(conditions),8);
% max_tscore = zeros(length(conditions),8);
% for i = 1:size(Lpos,1)
% this_tmatrix = spm_summarise(this_file,struct('def','sphere', 'spec',8, 'xyz',Lpos(i,:)'));
% mean_tscore(1,i) = mean(nonzeros(this_tmatrix));
% max_tscore(1,i) = max(nonzeros(this_tmatrix));
% end
% 
% 
% for j = 2:length(conditions)
%     this_file = ['/imaging/mlr/users/tc02/Holly_MMN/ICA_Denoise/VBM_MMN_amplitude_left A1_' conditions{j} '/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0001.nii'];
%     for i = 1:size(Lpos,1)
%         this_tmatrix = spm_summarise(this_file,struct('def','sphere', 'spec',8, 'xyz',Lpos(i,:)'));
%         mean_tscore(j,i) = mean(nonzeros(this_tmatrix));
%         max_tscore(j,i) = max(nonzeros(this_tmatrix));
%     end
% 
%     
%     %this_file = ['/imaging/mlr/users/tc02/Holly_MMN/ICA_Denoise/VBM_MMN_latency_left A1_' conditions{i} '/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0001.nii'];
% 
% end
% 
% 

% 
% jp_spm8_surfacerender2_version_tc(['/imaging/mlr/users/tc02/Holly_MMN/ICA_Denoise/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0001.nii'],'jet',cfg)
% savepath = ['./VBM/Con-PCA_001unc'];
% eval(['export_fig ' savepath '.png -transparent -m2.5'])
% close all
% 
% jp_spm8_surfacerender2_version_tc(['/imaging/mlr/users/tc02/Holly_MMN/ICA_Denoise/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0002.nii'],'jet',cfg)
% savepath = ['./VBM/Con-bvFTD_001unc'];
% eval(['export_fig ' savepath '.png -transparent -m2.5'])
% close all
% 
% jp_spm8_surfacerender2_version_tc(['/imaging/mlr/users/tc02/Holly_MMN/ICA_Denoise/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0003.nii'],'jet',cfg)
% savepath = ['./VBM/Con-nfvPPA_001unc'];
% eval(['export_fig ' savepath '.png -transparent -m2.5'])
% close all
% 
% jp_spm8_surfacerender2_version_tc(['/imaging/mlr/users/tc02/Holly_MMN/ICA_Denoise/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0007.nii'],'jet',cfg)
% savepath = ['./VBM/PCA-bvFTD_001unc'];
% eval(['export_fig ' savepath '.png -transparent -m2.5'])
% close all
% 
% jp_spm8_surfacerender2_version_tc(['/imaging/mlr/users/tc02/Holly_MMN/ICA_Denoise/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0009.nii'],'jet',cfg)
% savepath = ['./VBM/PCA-nfvPPA_001unc'];
% eval(['export_fig ' savepath '.png -transparent -m2.5'])
% close all
% 
% jp_spm8_surfacerender2_version_tc(['/imaging/mlr/users/tc02/Holly_MMN/ICA_Denoise/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0011.nii'],'jet',cfg)
% savepath = ['./VBM/bvFTD-nfvPPA_001unc'];
% eval(['export_fig ' savepath '.png -transparent -m2.5'])
% close all
% 
% jp_spm8_surfacerender2_version_tc(['/imaging/mlr/users/tc02/vespa/scans/PNFA_VBM/tom/Stats/ERF_factorial_full_group_vbm_TIVnormalised_agecovaried_smoothedmask/spmT_0001.nii'],'jet',cfg)
% savepath = ['./VBM/Con-nfvPPA_VESPA_001unc'];
% eval(['export_fig ' savepath '.png -transparent -m2.5'])
% close all
% 
% jp_spm8_surfacerender2_version_tc(['/imaging/mlr/users/tc02/vespa/scans/PNFA_VBM/tom/Stats/ERF_factorial_full_group_vbm_TIVnormalised_agecovaried_smoothedmask/spmT_0002.nii'],'jet',cfg)
% savepath = ['./VBM/Con-nfvPPA_VESPA_001unc_rev'];
% eval(['export_fig ' savepath '.png -transparent -m2.5'])
% close all
% 
% jp_spm8_surfacerender2_version_tc(['/imaging/mlr/users/tc02/Holly_MMN/ICA_Denoise/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0013.nii'],'jet',cfg)
% savepath = ['./VBM/Con-ADMCI_001unc'];
% eval(['export_fig ' savepath '.png -transparent -m2.5'])
% close all
% 
% jp_spm8_surfacerender2_version_tc(['/imaging/mlr/users/tc02/Holly_MMN/ICA_Denoise/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0015.nii'],'jet',cfg)
% savepath = ['./VBM/PCA-ADMCI_001unc'];
% eval(['export_fig ' savepath '.png -transparent -m2.5'])
% close all
% 
% jp_spm8_surfacerender2_version_tc(['/imaging/mlr/users/tc02/Holly_MMN/ICA_Denoise/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0017.nii'],'jet',cfg)
% savepath = ['./VBM/bvFTD-ADMCI_001unc'];
% eval(['export_fig ' savepath '.png -transparent -m2.5'])
% close all
% 
% jp_spm8_surfacerender2_version_tc(['/imaging/mlr/users/tc02/Holly_MMN/ICA_Denoise/VBM/VBM_stats/factorial_full_group_vbm_TIVnormalised_agecovaried_unsmoothedmask/spmT_0019.nii'],'jet',cfg)
% savepath = ['./VBM/nfvPPA-ADMCI_001unc'];
% eval(['export_fig ' savepath '.png -transparent -m2.5'])
% close all
