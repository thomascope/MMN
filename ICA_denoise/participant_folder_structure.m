% paths to different data
mf_path = '/imaging/rowe/archive/users/hp02/pnfa_mmn/maxfilter/';
mri_path = '/imaging/rowe/archive/users/hp02/pnfa_mmn/mri_scans/';
preproc_path = '/imaging/rowe/archive/users/hp02/pnfa_mmn/preprocessed/dvts_sep/';
LFP_path = 'imaging/tc02/scratch/LFPs_TA/';
dcm_path = '/imaging/rowe/archive/users/hp02/pnfa_mmn/dcm/CMC_DCM_all_subjs_together_camcanHCs/2019_TC_LFPs_customPriors_32models';

i=1;
% Participant details for patient MMN data analysis

Participant{i}.name = 'bvftd1_meg08_0379';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'sbvftd1_meg08_0379';
Participant{i}.MF = 'ftd_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'bvftd2_meg14_0046_cc520127';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'sbvftd2_meg14_0046_cc520127';
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'bvftd3_meg13_0363_cc510629';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'sbvftd3_meg13_0363_cc510629';
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'bvftd4_meg13_0479_cc420383';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI =  'sbvftd4_meg13_0479_cc420383';
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'bvftd5_meg13_0411_cc410084';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'sbvftd5_meg13_0411_cc410084'; 
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'bvftd6_meg14_0094';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'sbvftd6_meg14_0094'; 
Participant{i}.MF = 'ftd_0094_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'bvftd7_meg13_0411_cc510648';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'sbvftd7_meg13_0411_cc510648'; 
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'bvftd8_meg14_0051_cc420100';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'sbvftd8_meg14_0051_cc420100'; 
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'bvftd9_meg14_0067_cc420157';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'sbvftd9_meg14_0067_cc420157'; 
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'bvftd10_meg14_0114_vc4';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'sbvftd10_meg14_0114_vc4'; 
Participant{i}.MF = 'vespa_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'bvftd11_meg13_0506_cc410179';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'sbvftd11_meg13_0506_cc410179'; 
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'bvftd12_meg13_0403_cc510639';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'sbvftd12_meg13_0403_cc510639'; 
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'bvftd13_meg13_0393_cc520065';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'sbvftd13_meg13_0393_cc520065'; 
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'bvftd14_meg13_0480_cc420261';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'sbvftd14_meg13_0480_cc420261'; 
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'bvftd15_meg13_0338_cc412021';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'sbvftd15_meg13_0338_cc412021'; 
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'bvftd16_meg08_0368';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'sbvftd16_meg08_0368'; 
Participant{i}.MF = 'ftd_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'bvftd17_meg12_0040';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'sbvftd17_meg12_0040'; 
Participant{i}.MF = 'FTD_40_MMN_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'bvftd18_meg13_0317_cc520745';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'sbvftd18_meg13_0317_cc520745'; 
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'bvftd19_meg14_0060_cc520253';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'sbvftd19_meg14_0060_cc520253'; 
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'bvftd20_meg13_0401_cc420198';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'sbvftd20_meg13_0401_cc420198'; 
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'bvftd21_meg14_0096_vc3';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'sbvftd21_meg14_0096_vc3'; 
Participant{i}.MF = 'vespa_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'bvftd22_meg14_0271_cc520097';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'sbvftd22_meg14_0271_cc520097'; 
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'bvftd23_meg14_0171_cc610496';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'sbvftd23_meg14_0171_cc610496'; 
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pca1_meg14_0044_cc510355';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'spca1_meg14_0044_cc510355'; 
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pca2_meg09_0087';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'spca2_meg09_0087'; 
Participant{i}.MF = 'ftd_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pca3_meg08_0365';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'spca3_meg08_0365'; 
Participant{i}.MF = 'ftd_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pca4_meg12_0520';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'spca4_meg12_521'; 
Participant{i}.MF = 'ftd_mmn_520_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pca5_meg12_0036';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'spca5_meg12_0036'; 
Participant{i}.MF = 'FTD_36_MMN_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pca6_meg14_0056_cc410284';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'spca6_meg14_0056_cc410284'; 
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pca7_meg14_0009_cc420435';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'spca7_meg14_0009_cc420435'; 
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pca8_meg14_0249_cc510259';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'spca8_meg14_0249_cc510259'; 
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pca9_meg14_0048_cc512003';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'spca9_meg14_0048_cc512003'; 
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pca10_meg11_0267';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'spca10_meg11_0267'; 
Participant{i}.MF = 'ftd_mmn_0267_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pca11_meg13_0372_cc710446';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'spca11_meg13_0372_cc710446'; 
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pca12_meg14_0235_cc420486';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'spca12_meg14_0235_cc420486'; 
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pca13_meg13_0541_cc410390';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'spca13_meg13_0541_cc410390'; 
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pca14_meg12_0031';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'spca14_meg12_0031'; 
Participant{i}.MF = 'FTD_31_MMN_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pca15_meg14_0016_cc610405';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI = 'spca15_meg14_0016_cc610405'; 
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pnfa1_meg14_0215_cc620572';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI =  'spnfa1_meg14_0215_cc620572';
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pnfa2_meg14_0011_cc420162';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI =  'spnfa2_meg14_0011_cc420162';
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pnfa3_meg14_0086_cc520055';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI =  'spnfa3_meg14_0086_cc520055';
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pnfa4_meg13_0404_cc620592';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI =  'spnfa4_meg13_0404_cc620592';
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pnfa5_meg08_0273';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI =  'spnfa5_meg08_0273';
Participant{i}.MF = 'ftdmmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pnfa6_meg14_0120_cc610061';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI =  'spnfa6_meg14_0120_cc610061';
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pnfa7_meg13_0418_cc610288';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI =  'spnfa7_meg13_0418_cc610288';
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pnfa8_meg14_0428_cc721291';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI =  'spnfa8_meg14_0428_cc721291';
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pnfa9_meg14_0441_cc610932';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI =  'spnfa9_meg14_0441_cc610932';
Participant{i}.MF = 'mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'pnfa10_meg08_0378';
Participant{i}.groupfolder = 'matched_HCs';
Participant{i}.diag = 'Control';
Participant{i}.MRI =  'spnfa10_meg08_0378';
Participant{i}.MF = 'ftd_mmn_raw_trans1stdef.fif';
i=i+1;


% bvFTDs

Participant{i}.name = 'meg11_0179';%_nomri
Participant{i}.groupfolder = 'bvftd';
Participant{i}.diag = 'bvFTD';
Participant{i}.MRI = 'single_subj_T1';
Participant{i}.MF = '110708/ftd_mmn_0179_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg11_0238_2';
Participant{i}.groupfolder = 'bvftd';
Participant{i}.diag = 'bvFTD';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '110915/ftd_mmn_2_0238_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg11_0249';
Participant{i}.groupfolder = 'bvftd';
Participant{i}.diag = 'bvFTD';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '110927/ftd_mmn_1_2_3_0249_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg11_0270';
Participant{i}.groupfolder = 'bvftd';
Participant{i}.diag = 'bvFTD';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '111018/ftd_mmn_0270_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg12_0060';
Participant{i}.groupfolder = 'bvftd';
Participant{i}.diag = 'bvFTD';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '120307/ftd_0060_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg12_0072';
Participant{i}.groupfolder = 'bvftd';
Participant{i}.diag = 'bvFTD';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '120314/ftd_0072_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg12_0092';
Participant{i}.groupfolder = 'bvftd';
Participant{i}.diag = 'bvFTD';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '120323/ftd_0092_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg12_0161';%_nomri
Participant{i}.groupfolder = 'bvftd';
Participant{i}.diag = 'bvFTD';
Participant{i}.MRI = 'single_subj_T1';
Participant{i}.MF = '120426/ftd_mmn_0161_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg12_0228';%_nomri
Participant{i}.groupfolder = 'bvftd';
Participant{i}.diag = 'bvFTD';
Participant{i}.MRI = 'single_subj_T1';
Participant{i}.MF = '120523/ftd_0228_mmn2_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg12_0314';%_no_mri
Participant{i}.groupfolder = 'bvftd';
Participant{i}.diag = 'bvFTD';
Participant{i}.MRI = 'single_subj_T1';
Participant{i}.MF = '120709/ftd_0314_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg12_0366_2';
Participant{i}.groupfolder = 'bvftd';
Participant{i}.diag = 'bvFTD';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '120723/ftd_mmn3_0366_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg12_0389';
Participant{i}.groupfolder = 'bvftd';
Participant{i}.diag = 'bvFTD';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '120808/ftd_mmn_0389_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg12_0496_2';
Participant{i}.groupfolder = 'bvftd';
Participant{i}.diag = 'bvFTD';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '121019/ftd_mmn_1_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg12_0504';
Participant{i}.groupfolder = 'bvftd';
Participant{i}.diag = 'bvFTD';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '121025/ftd_mmn_0504_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg13_0016';%_nomri
Participant{i}.groupfolder = 'bvftd';
Participant{i}.diag = 'bvFTD';
Participant{i}.MRI = 'single_subj_T1';
Participant{i}.MF = '130125/ftd_mmn_1_2_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg13_0162';
Participant{i}.groupfolder = 'bvftd';
Participant{i}.diag = 'bvFTD';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '130405/ftd_mmn_0162_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg13_0279';
Participant{i}.groupfolder = 'bvftd';
Participant{i}.diag = 'bvFTD';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '130605/ftd_0279_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg13_0410';
Participant{i}.groupfolder = 'bvftd';
Participant{i}.diag = 'bvFTD';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '130919/ftd_13_0410_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg13_0454';
Participant{i}.groupfolder = 'bvftd';
Participant{i}.diag = 'bvFTD';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '131017/ftd_0454_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg14_0143';
Participant{i}.groupfolder = 'bvftd';
Participant{i}.diag = 'bvFTD';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '140331/ftd_0143_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg14_0315';
Participant{i}.groupfolder = 'bvftd';
Participant{i}.diag = 'bvFTD';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '140724/ftd_14_0315_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg14_0559';
Participant{i}.groupfolder = 'bvftd';
Participant{i}.diag = 'bvFTD';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '141215/140559_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg15_0061';%_nomri
Participant{i}.groupfolder = 'bvftd';
Participant{i}.diag = 'bvFTD';
Participant{i}.MRI = 'single_subj_T1';
Participant{i}.MF = '150319/ftd_15_0061_mmn_raw_trans1stdef.fif';
i=i+1;

% PCA
Participant{i}.name = 'meg13_0220';
Participant{i}.groupfolder = 'pca';
Participant{i}.diag = 'pca';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '130501/ftd_mmn_0220_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg13_0225';
Participant{i}.groupfolder = 'pca';
Participant{i}.diag = 'pca';
Participant{i}.MRI = '22032_structural';
Participant{i}.MF = '130502/ftd_0225_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg13_0236';
Participant{i}.groupfolder = 'pca';
Participant{i}.diag = 'pca';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '130508/ftd_0236_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg13_0277';
Participant{i}.groupfolder = 'pca';
Participant{i}.diag = 'pca';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '130603/ftd_0277_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg13_0284';
Participant{i}.groupfolder = 'pca';
Participant{i}.diag = 'pca';
Participant{i}.MRI = '24192_0002';
Participant{i}.MF = '130607/ftd_0284_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg13_0300';
Participant{i}.groupfolder = 'pca';
Participant{i}.diag = 'pca';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '130628/ftd_0300_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg13_0303';
Participant{i}.groupfolder = 'pca';
Participant{i}.diag = 'pca';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '130702/ftd_303_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg13_0324';
Participant{i}.groupfolder = 'pca';
Participant{i}.diag = 'pca';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '130723/ftd_0324_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg13_0356';
Participant{i}.groupfolder = 'pca';
Participant{i}.diag = 'pca';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '130814/ftd_0356_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg13_0437';
Participant{i}.groupfolder = 'pca';
Participant{i}.diag = 'pca';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '131008/ftd_0437_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg14_0061';%_nomri
Participant{i}.groupfolder = 'pca';
Participant{i}.diag = 'pca';
Participant{i}.MRI = 'single_subj_T1';
Participant{i}.MF = '140217/ftd_0061_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg14_0199';%_nomri
Participant{i}.groupfolder = 'pca';
Participant{i}.diag = 'pca';
Participant{i}.MRI = 'single_subj_T1';
Participant{i}.MF = '140430/ftd_0199_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg14_0287';
Participant{i}.groupfolder = 'pca';
Participant{i}.diag = 'pca';
Participant{i}.MRI = '24599_0002';
Participant{i}.MF = '140624/ftd_0287_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg14_0327';
Participant{i}.groupfolder = 'pca';
Participant{i}.diag = 'pca';
Participant{i}.MRI = 'structural';
Participant{i}.MF = '140804/0327_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg14_0333';%_nomri
Participant{i}.groupfolder = 'pca';
Participant{i}.diag = 'pca';
Participant{i}.MRI = 'single_subj_T1';
Participant{i}.MF = '140808/ftd_meg14_0333_MMN_raw_trans1stdef.fif';
i=i+1;

% nfvPPA


Participant{i}.name = 'meg14_0085_vp1';
Participant{i}.groupfolder = 'vespa';
Participant{i}.diag = 'nfvppa';
Participant{i}.MRI = 'vp1_Structural';
Participant{i}.MF = '140228/vespa_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg14_0087_vp2';
Participant{i}.groupfolder = 'vespa';
Participant{i}.diag = 'nfvppa';
Participant{i}.MRI = 'vp2_Structural';
Participant{i}.MF = '140303/vespa_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg14_0121_vp5';
Participant{i}.groupfolder = 'vespa';
Participant{i}.diag = 'nfvppa';
Participant{i}.MRI = 'vp5_Structural';
Participant{i}.MF = '140318/vespa_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg14_0130_vp6';
Participant{i}.groupfolder = 'vespa';
Participant{i}.diag = 'nfvppa';
Participant{i}.MRI = 'vp6_Structural';
Participant{i}.MF = '140324/vespa_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg14_0135_vp7';
Participant{i}.groupfolder = 'vespa';
Participant{i}.diag = 'nfvppa';
Participant{i}.MRI = 'avg152T1';
Participant{i}.MF = '140325/vespa_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg14_0150_vp8';
Participant{i}.groupfolder = 'vespa';
Participant{i}.diag = 'nfvppa';
Participant{i}.MRI = 'vp8_Structural';
Participant{i}.MF = '140403/vespa_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg14_0184_vp9';
Participant{i}.groupfolder = 'vespa';
Participant{i}.diag = 'nfvppa';
Participant{i}.MRI = 'vp9_Structural';
Participant{i}.MF = '140424/vespa_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg14_0205_vp10';
Participant{i}.groupfolder = 'vespa';
Participant{i}.diag = 'nfvppa';
Participant{i}.MRI = 'vp10_Structural';
Participant{i}.MF = '140502/vespa_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg14_0222_vp11';
Participant{i}.groupfolder = 'vespa';
Participant{i}.diag = 'nfvppa';
Participant{i}.MRI = 'vp11_Structural';
Participant{i}.MF = '140513/vespa_mmn_raw_trans1stdef.fif';
i=i+1;

Participant{i}.name = 'meg14_0506_vp12';
Participant{i}.groupfolder = 'vespa';
Participant{i}.diag = 'nfvppa';
Participant{i}.MRI = 'vp12_Structural';
Participant{i}.MF = '140506/vespa_mmn_raw_trans1stdef.fif';
i=i+1;

% %MCIs
% addpath('/imaging/rowe/archive/users/hp02/pnfa_mmn/preprocessed/For_Thomas_dvts_sep/mri_scans');
% copy_mci_mris;% This script was used to copy mris, but it also contains mri info useful to quickly make MCI struct
% for pp = 1:length(mci_pats)
%     Participant{i}.name = mci_pats{pp,1};
%     Participant{i}.groupfolder = 'MCI';
%     Participant{i}.diag = 'MCI';
%     
%     [mri_path, mri_name, mri_ext] = fileparts(pp_mri_path{pp,1})
%     Participant{i}.MRI = mri_name;
%     %Participant{i}.MF = not sure what these names are yet
%     i = i+1;
% end
% 
% %% for folders
% 
% for pp = 1:length(Participant)
%     
% end