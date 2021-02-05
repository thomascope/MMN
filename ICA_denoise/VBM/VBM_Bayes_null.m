VBM_path = '/imaging/tc02/Holly_MMN/ICA_Denoise/VBM/VBM_stats/Bayesian';

%Contrasts:
this_fname = 'Con-PCA_null';
this_bfname = 'logbf_0021';
spm_bms_test_null_nameoutput([VBM_path '/' this_bfname '.nii'], this_fname)

this_fname = 'Con-bvFTD_null';
this_bfname = 'logbf_0022';
spm_bms_test_null_nameoutput([VBM_path '/' this_bfname '.nii'], this_fname)

this_fname = 'Con-nfvPPA_null';
this_bfname = 'logbf_0023';
spm_bms_test_null_nameoutput([VBM_path '/' this_bfname '.nii'], this_fname)

this_fname = 'Con-ADMCI_null';
this_bfname = 'logbf_0024';
spm_bms_test_null_nameoutput([VBM_path '/' this_bfname '.nii'], this_fname)

