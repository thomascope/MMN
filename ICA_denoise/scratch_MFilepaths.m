for i = 1:length(Participant)
    this_MF_source{i} = {};
    if iscell(Participant{i}.MF)
        for j = 1:length(Participant{i}.MF)
            if exist(Participant{i}.MF{j},'file')
                this_MF_source{i}{end+1} = Participant{i}.MF{j};
            else
                warning(['No Maxfilter source found for Participant ' num2str(i)])
            end
        end
    elseif exist(Participant{i}.MF,'file')
        this_MF_source{i}{1} = Participant{i}.MF
    elseif exist(['/imaging/rowe/archive/users/hp02/pnfa_mmn/maxfilter/' Participant{i}.groupfolder '/' Participant{i}.name '/' Participant{i}.MF],'file')
        this_MF_source{i}{1} = ['/imaging/rowe/archive/users/hp02/pnfa_mmn/maxfilter/' Participant{i}.groupfolder '/' Participant{i}.name '/' Participant{i}.MF];
    elseif exist(['/imaging/rowe/archive/users/hp02/pnfa_mmn/maxfilter/' Participant{i}.groupfolder '/' Participant{i}.name(1:end-2) '/' Participant{i}.MF],'file')
        this_MF_source{i}{1} = ['/imaging/rowe/archive/users/hp02/pnfa_mmn/maxfilter/' Participant{i}.groupfolder '/' Participant{i}.name(1:end-2) '/' Participant{i}.MF];
    else
        this_subfolder = ls(['/imaging/rowe/archive/users/hp02/pnfa_mmn/maxfilter/' Participant{i}.groupfolder '/' Participant{i}.name '/']);
        if exist(['/imaging/rowe/archive/users/hp02/pnfa_mmn/maxfilter/' Participant{i}.groupfolder '/' Participant{i}.name '/' this_subfolder '/' Participant{i}.MF],'file')
            this_MF_source{i}{1} = ['/imaging/rowe/archive/users/hp02/pnfa_mmn/maxfilter/' Participant{i}.groupfolder '/' Participant{i}.name '/' this_subfolder '/' Participant{i}.MF];
        else
            warning(['No Maxfilter source found for Participant ' num2str(i)])
        end
    end
end


%% fif files
%bvFTD
ftd_fif_path = {{'/megdata/cbu/ftd/meg11_0179/110708/ftd_mmn_0179_raw.fif'};
{'/megdata/cbu/ftd/meg11_0238/110915/ftd_mmn_1_0238_raw.fif','/megdata/cbu/ftd/meg11_0238/110915/ftd_mmn_2_0238_raw.fif'};
{'/megdata/cbu/ftd/meg11_0249/110927/ftd_mmn_1_2_3_0249_raw.fif'};
{'/megdata/cbu/ftd/meg11_0270/111018/ftd_mmn_0270_raw.fif'};
{'/megdata/cbu/ftd/meg12_0060/120307/ftd_0060_mmn_raw.fif'};
{'/megdata/cbu/ftd/meg12_0072/120314/ftd_0072_mmn_raw.fif'};
{'/megdata/cbu/ftd/meg12_0092/120323/ftd_0092_mmn_raw.fif'};
{'/megdata/cbu/ftd/meg12_0161/120426/ftd_mmn_0161_raw.fif'};
{'/megdata/cbu/ftd/meg12_0228/120523/ftd_0228_mmn1_raw.fif','/megdata/cbu/ftd/meg12_0228/120523/ftd_0228_mmn2_raw.fif'};
{'/megdata/cbu/ftd/meg12_0314/120709/ftd_0314_mmn_raw.fif'};
{'/megdata/cbu/ftd/meg12_0366/120723/ftd_mmn_0366_raw.fif','/megdata/cbu/ftd/meg12_0366/120723/ftd_mmn3_0366_raw.fif'};
{'/megdata/cbu/ftd/meg12_0389/120808/ftd_mmn_0389_raw.fif'};
{'/megdata/cbu/ftd/meg12_0496/121019/ftd_mmn_first_try_raw.fif','/megdata/cbu/ftd/meg12_0496/121019/ftd_mmn_1_raw.fif'};
{'/megdata/cbu/ftd/meg12_0504/121025/ftd_mmn_0504_raw.fif'};
%{'/megdata/cbu/ftd/meg/ftd_mmn_0519_raw.fif'};
{'/megdata/cbu/ftd/meg13_0016/130125/ftd_mmn_1_2_raw.fif'};
{'/megdata/cbu/ftd/meg13_0162/130405/ftd_mmn_0162_raw.fif'};
{'/megdata/cbu/ftd/meg13_0279/130605/ftd_0279_mmn_raw.fif'};
%{'/megdata/cbu/ftd/meg/ftd_mmn_0315_raw.fif'};
{'/megdata/cbu/ftd/meg13_0410/130919/ftd_13_0410_mmn_raw.fif'};
{'/megdata/cbu/ftd/meg13_0454/131017/ftd_0454_mmn_raw.fif'};
{'/megdata/cbu/ftd/meg14_0143/140331/ftd_0143_mmn_raw.fif'};
{'/megdata/cbu/ftd/meg14_0315/140724/ftd_14_0315_mmn_raw.fif'}; 
{'/megdata/cbu/ftd/meg14_0559/141215/140559_mmn_raw.fif'};
{'/megdata/cbu/ftd/meg15_0061/150319/ftd_15_0061_mmn_raw.fif'}}; 
ftd_pats = {{'ftd_mmn_0179_raw.fif'};
{'ftd_mmn_1_0238_raw.fif','ftd_mmn_2_0238_raw.fif'};
{'ftd_mmn_1_2_3_0249_raw.fif'};
{'ftd_mmn_0270_raw.fif'};
{'ftd_0060_mmn_raw.fif'};
{'ftd_0072_mmn_raw.fif'};
{'ftd_0092_mmn_raw.fif'};
{'ftd_mmn_0161_raw.fif'};
{'ftd_0228_mmn1_raw.fif','ftd_0228_mmn2_raw.fif'};
{'ftd_0314_mmn_raw.fif'};
{'ftd_mmn_0366_raw.fif','ftd_mmn3_0366_raw.fif'};
{'ftd_mmn_0389_raw.fif'};
{'ftd_mmn_first_try_raw.fif','ftd_mmn_1_raw.fif'};
{'ftd_mmn_0504_raw.fif'};
%{'ftd_mmn_0519_raw.fif'};
{'ftd_mmn_1_2_raw.fif'};
{'ftd_mmn_0162_raw.fif'};
{'ftd_0279_mmn_raw.fif'};
{'ftd_13_0410_mmn_raw.fif'};
{'ftd_mmn_0315_raw.fif'};
{'ftd_0454_mmn_raw.fif'};
{'ftd_0143_mmn_raw.fif'};
{'ftd_14_0315_mmn_raw.fif'}; 
{'140559_mmn_raw.fif'};
{'ftd_15_0061_mmn_raw.fif'}}; 

% PCA
pca_pats = [{'ftd_mmn_0220_raw.fif'};
{'ftd_0225_mmn_raw.fif'};
{'ftd_0236_mmn_raw.fif'}; 
{'ftd_0277_mmn_raw.fif'};
{'ftd_0284_mmn_raw.fif'};
{'ftd_0300_mmn_raw.fif'};
{'ftd_303_mmn_raw.fif'};
{'ftd_0324_mmn_raw.fif'}; 
{'ftd_0356_mmn_raw.fif'}; 
{'ftd_0437_mmn_raw.fif'}; 
%{'ftd_0526_mmn_raw.fif'};
{'ftd_0061_mmn_raw.fif'};
{'ftd_0199_mmn_raw.fif'};
{'ftd_0287_mmn_raw.fif'}; 
{'0327_mmn_raw.fif'}; 
{'ftd_meg14_0333_MMN_raw.fif'}]; 

pca_fif_path = [{'/megdata/cbu/ftd/meg13_0220/130501/ftd_mmn_0220_raw.fif'};
{'/megdata/cbu/ftd/meg13_0225/130502/ftd_0225_mmn_raw.fif'};
{'/megdata/cbu/ftd/meg13_0236/130508/ftd_0236_mmn_raw.fif'}; 
{'/megdata/cbu/ftd/meg13_0277/130603/ftd_0277_mmn_raw.fif'};
{'/megdata/cbu/ftd/meg13_0284/130607/ftd_0284_mmn_raw.fif'};
{'/megdata/cbu/ftd/meg13_0300/130628/ftd_0300_mmn_raw.fif'};
{'/megdata/cbu/ftd/meg13_0303/130702/ftd_303_mmn_raw.fif'};
{'/megdata/cbu/ftd/meg13_0324/130723/ftd_0324_mmn_raw.fif'}; 
{'/megdata/cbu/ftd/meg13_0356/130814/ftd_0356_mmn_raw.fif'}; 
{'/megdata/cbu/ftd/meg13_0437/131008/ftd_0437_mmn_raw.fif'}; 
%{'/megdata/cbu/ftd/ftd_0526_mmn_raw.fif'};
{'/megdata/cbu/ftd/meg14_0061/140217/ftd_0061_mmn_raw.fif'};
{'/megdata/cbu/ftd/meg14_0199/140430/ftd_0199_mmn_raw.fif'};
{'/megdata/cbu/ftd/meg14_0287/140624/ftd_0287_mmn_raw.fif'}; 
{'/megdata/cbu/ftd/meg14_0327/140804/0327_mmn_raw.fif'}; 
{'/megdata/cbu/ftd/meg14_0333/140808/ftd_meg14_0333_MMN_raw.fif'}]; 

% VESPA - from Thomas
vespa_fif_path = [{'/megdata/cbu/vespa/meg14_0085_vp1/140228/vespa_mmn_raw.fif';
    '/megdata/cbu/vespa/meg14_0087_vp2/140303/vespa_mmn_raw.fif';
    '/megdata/cbu/vespa/meg14_0121_vp5/140318/vespa_mmn_raw.fif';
    '/megdata/cbu/vespa/meg14_0130_vp6/140324/vespa_mmn_raw.fif';
    '/megdata/cbu/vespa/meg14_0135_vp7/140325/vespa_mmn_raw.fif';
    '/megdata/cbu/vespa/meg14_0150_vp8/140403/vespa_mmn_raw.fif';
    '/megdata/cbu/vespa/meg14_0184_vp9/140424/vespa_mmn_raw.fif';
    '/megdata/cbu/vespa/meg14_0205_vp10/140502/vespa_mmn_raw.fif';
    '/megdata/cbu/vespa/meg14_0222_vp11/140513/vespa_mmn_raw.fif';
    '/imaging/mlr/users/tc02/vespa/preprocess/meg14_0506_vp12/MMN+Rest/mmn_raw_sss.fif'}];

match_HCs = {
    % bvFTD matched controls:
    '/megdata/cbu/ftd/meg08_0379/080829/'
    '/megdata/camcan/camcan_two/meg14_0046_cc520127/140206/'
    '/megdata/camcan/camcan_two/meg13_0363_cc510629/130820/'
    '/megdata/camcan/camcan_two/meg13_0479_cc420383/131104/'
    '/megdata/camcan/camcan_two/meg13_0411_cc410084/130920/'
    '/megdata/cbu/ftd/meg14_0094/140306/'
    '/megdata/camcan/camcan_two/meg13_0411_cc510648/130919/'
    '/megdata/camcan/camcan_two/meg14_0051_cc420100/140210/'
    '/megdata/camcan/camcan_two/meg14_0067_cc420157/140220/'
    '/megdata/cbu/vespa/meg14_0114_vc4/140314/'
    '/megdata/camcan/camcan_two/meg13_0506_cc410179/131123/'
    '/megdata/camcan/camcan_two/meg13_0403_cc510639/130913/'
    '/megdata/camcan/camcan_two/meg13_0393_cc520065/130907/'
    '/megdata/camcan/camcan_two/meg13_0480_cc420261/131105/'
    '/megdata/camcan/camcan_two/meg13_0338_cc412021/130802/'
    '/megdata/cbu/ftd/meg08_0368/080822/'
    '/megdata/cbu/ftd/meg12_0040/120224/'
    '/megdata/camcan/camcan_two/meg13_0317_cc520745/130718/'
    '/megdata/camcan/camcan_two/meg14_0060_cc520253/140217/'
    '/megdata/camcan/camcan_two/meg13_0401_cc420198/130912/'
    '/megdata/cbu/vespa/meg14_0096_vc3/140307/'
    '/megdata/camcan/camcan_two/meg14_0271_cc520097/140610/'
    '/megdata/camcan/camcan_two/meg14_0171_cc610496/140415/'
    
    %PCA matched controls:
    '/megdata/camcan/camcan_two/meg14_0044_cc510355/140204/'
    '/megdata/cbu/ftd/meg09_0087/090326/'
    '/megdata/cbu/ftd/meg08_0365/080821/'
    '/megdata/cbu/ftd/meg12_0520/121106/'
    '/megdata/cbu/ftd/meg12_0036/120223/'
    '/megdata/camcan/camcan_two/meg14_0056_cc410284/140213/'
    '/megdata/camcan/camcan_two/meg14_0009_cc420435/140117/'
    '/megdata/camcan/camcan_two/meg14_0249_cc510259/140528/'
    '/megdata/camcan/camcan_two/meg14_0048_cc512003/140207/'
    '/megdata/cbu/ftd/meg11_0267/111014/'
    '/megdata/camcan/camcan_two/meg13_0372_cc710446/130823/'
    '/megdata/camcan/camcan_two/meg14_0235_cc420486/140519/'
    '/megdata/camcan/camcan_two/meg13_0541_cc410390/131216/'
    '/megdata/cbu/ftd/meg12_0031/120221/'
    '/megdata/camcan/camcan_two/meg14_0016_cc610405/140121/'
    
    % PNFA matched controls
    '/megdata/camcan/camcan_two/meg14_0215_cc620572/140509/'
    '/megdata/camcan/camcan_two/meg14_0011_cc420162/140117/'
    '/megdata/camcan/camcan_two/meg14_0086_cc520055/140303/'
    '/megdata/camcan/camcan_two/meg13_0404_cc620592/130916/'
    '/megdata/cbu/ftd/meg08_0273/080609/'
    '/megdata/camcan/camcan_two/meg14_0120_cc610061/140318/'
    '/megdata/camcan/camcan_two/meg13_0418_cc610288/130926/'
    '/megdata/camcan/camcan_two/meg14_0428_cc721291/140127/'
    '/megdata/camcan/camcan_two/meg14_0441_cc610932/140516/'
    '/megdata/cbu/ftd/meg08_0378/080829/'
    };

match_HCs_names = {
    % bvFTD matched controls:
    'bvftd1_meg08_0379';
    'bvftd2_meg14_0046_cc520127';
    'bvftd3_meg13_0363_cc510629';
    'bvftd4_meg13_0479_cc420383';
    'bvftd5_meg13_0411_cc410084';
    'bvftd6_meg14_0094';
    'bvftd7_meg13_0411_cc510648';
    'bvftd8_meg14_0051_cc420100';
    'bvftd9_meg14_0067_cc420157';
    'bvftd10_meg14_0114_vc4';
    'bvftd11_meg13_0506_cc410179';
    'bvftd12_meg13_0403_cc510639';
    'bvftd13_meg13_0393_cc520065';
    'bvftd14_meg13_0480_cc420261';
    'bvftd15_meg13_0338_cc412021';
    'bvftd16_meg08_0368';
    'bvftd17_meg12_0040';
    'bvftd18_meg13_0317_cc520745';
    'bvftd19_meg14_0060_cc520253';
    'bvftd20_meg13_0401_cc420198';
    'bvftd21_meg14_0096_vc3';
    'bvftd22_meg14_0271_cc520097';
    'bvftd23_meg14_0171_cc610496';
    
    %PCA matched controls:
    'pca1_meg14_0044_cc510355';
    'pca2_meg09_0087';
    'pca3_meg08_0365';
    'pca4_meg12_0520';
    'pca5_meg12_0036';
    'pca6_meg14_0056_cc410284';
    'pca7_meg14_0009_cc420435';
    'pca8_meg14_0249_cc510259';
    'pca9_meg14_0048_cc512003';
    'pca10_meg11_0267';
    'pca11_meg13_0372_cc710446';
    'pca12_meg14_0235_cc420486';
    'pca13_meg13_0541_cc410390';
    'pca14_meg12_0031';
    'pca15_meg14_0016_cc610405';
    
    % PNFA matched controls
    'pnfa1_meg14_0215_cc620572';
    'pnfa2_meg14_0011_cc420162';
    'pnfa3_meg14_0086_cc520055';
    'pnfa4_meg13_0404_cc620592';
    'pnfa5_meg08_0273';
    'pnfa6_meg14_0120_cc610061';
    'pnfa7_meg13_0418_cc610288';
    'pnfa8_meg14_0428_cc721291';
    'pnfa9_meg14_0441_cc610932';
    'pnfa10_meg08_0378';
    };

% Order of 08/09 files (2) and 2010-19 files (import for dir function)
matched_HC_date_order = [2	1	1	1	1	1	1	1	1	1	1	1	1	1	1	2	1	1	1	1	1	1	1		1	2	2	1	1	1	1	1	1	1	1	1	1	1	1			1	1	1	1	2	1	1	1	1	2];

% Order of mmn fif file names:
mmn_fif_names_mHCs = {
% bvFTD
'ftd_mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'ftd_0094_mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'vespa_mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'ftd_mmn_raw_trans1stdef.fif';
'FTD_40_MMN_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'vespa_mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
% PCA
'mmn_raw_trans1stdef.fif';
'ftd_mmn_raw_trans1stdef.fif';
'ftd_mmn_raw_trans1stdef.fif';
'ftd_mmn_520_raw_trans1stdef.fif';
'FTD_36_MMN_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'ftd_mmn_0267_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'FTD_31_MMN_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';

% PNFA
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'ftdmmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'ftd_mmn_raw_trans1stdef.fif';
};


%% fif names after maxfilter:

% Order of mmn fif file names:
mmn_fif_names_mHCs_transdef = {
% bvFTD
'ftd_mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'ftd_0094_mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'vespa_mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'ftd_mmn_raw_trans1stdef.fif';
'FTD_40_MMN_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'vespa_mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
% PCA
'mmn_raw_trans1stdef.fif';
'ftd_mmn_raw_trans1stdef.fif';
'ftd_mmn_raw_trans1stdef.fif';
'ftd_mmn_520_raw_trans1stdef.fif';
'FTD_36_MMN_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'ftd_mmn_0267_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'FTD_31_MMN_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';

% PNFA
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'ftdmmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'mmn_raw_trans1stdef.fif';
'ftd_mmn_raw_trans1stdef.fif';
};

%% Preprocessed files
preproc_path = '/imaging/rowe/archive/users/hp02/pnfa_mmn/preprocessed/dvts_sep/';

preproc_matched_hcs = {
  
'matched_HCs/bvftd1_meg08_0379/fmraedfffbvftd1_meg08_0379.nii';
'matched_HCs/bvftd2_meg14_0046_cc520127/fmraedfffbvftd2_meg14_0046_cc520127.nii'
'matched_HCs/bvftd3_meg13_0363_cc510629/fmraedfffbvftd3_meg13_0363_cc510629.nii'
'matched_HCs/bvftd4_meg13_0479_cc420383/fmraedfffbvftd4_meg13_0479_cc420383.nii'
'matched_HCs/bvftd5_meg13_0411_cc410084/fmraedfffbvftd5_meg13_0411_cc410084.nii'
'matched_HCs/bvftd6_meg14_0094/fmraedfffbvftd6_meg14_0094.nii'
'matched_HCs/bvftd7_meg13_0411_cc510648/fmraedfffbvftd7_meg13_0411_cc510648.nii'
'matched_HCs/bvftd8_meg14_0051_cc420100/fmraedfffbvftd8_meg14_0051_cc420100.nii'
'matched_HCs/bvftd9_meg14_0067_cc420157/fmraedfffbvftd9_meg14_0067_cc420157.nii'
'matched_HCs/bvftd10_meg14_0114_vc4/fmraedfffbvftd10_meg14_0114_vc4.nii'
'matched_HCs/bvftd11_meg13_0506_cc410179/fmraedfffbvftd11_meg13_0506_cc410179.nii'
'matched_HCs/bvftd12_meg13_0403_cc510639/fmraedfffbvftd12_meg13_0403_cc510639.nii'
'matched_HCs/bvftd13_meg13_0393_cc520065/fmraedfffbvftd13_meg13_0393_cc520065.nii'
'matched_HCs/bvftd14_meg13_0480_cc420261/fmraedfffbvftd14_meg13_0480_cc420261.nii'
'matched_HCs/bvftd15_meg13_0338_cc412021/fmraedfffbvftd15_meg13_0338_cc412021.nii'
'matched_HCs/bvftd16_meg08_0368/fmraedfffbvftd16_meg08_0368.nii'
'matched_HCs/bvftd17_meg12_0040/fmraedfffbvftd17_meg12_0040.nii'
'matched_HCs/bvftd18_meg13_0317_cc520745/fmraedfffbvftd18_meg13_0317_cc520745.nii'
'matched_HCs/bvftd19_meg14_0060_cc520253/fmraedfffbvftd19_meg14_0060_cc520253.nii'
'matched_HCs/bvftd20_meg13_0401_cc420198/fmraedfffbvftd20_meg13_0401_cc420198.nii'
'matched_HCs/bvftd21_meg14_0096_vc3/fmraedfffbvftd21_meg14_0096_vc3.nii'
'matched_HCs/bvftd22_meg14_0271_cc520097/fmraedfffbvftd22_meg14_0271_cc520097.nii'
'matched_HCs/bvftd23_meg14_0171_cc610496/fmraedfffbvftd23_meg14_0171_cc610496.nii'

'matched_HCs/pca1_meg14_0044_cc510355/fmraedfffpca1_meg14_0044_cc510355.nii'
'matched_HCs/pca2_meg09_0087/fmraedfffpca2_meg09_0087.nii'
'matched_HCs/pca3_meg08_0365/fmraedfffpca3_meg08_0365.nii'
'matched_HCs/pca4_meg12_0520/fmraedfffpca4_meg12_0520.nii'
'matched_HCs/pca5_meg12_0036/fmraedfffpca5_meg12_0036.nii'
'matched_HCs/pca6_meg14_0056_cc410284/fmraedfffpca6_meg14_0056_cc410284/.nii'
'matched_HCs/pca7_meg14_0009_cc420435/fmraedfffpca7_meg14_0009_cc420435.nii'
'matched_HCs/pca8_meg14_0249_cc510259/fmraedfffpca8_meg14_0249_cc510259.nii'
'matched_HCs/pca9_meg14_0048_cc512003/fmraedfffpca9_meg14_0048_cc512003.nii'
'matched_HCs/pca10_meg11_0267/fmraedfffpca10_meg11_0267.nii'
'matched_HCs/pca11_meg13_0372_cc710446/fmraedfffpca11_meg13_0372_cc710446.nii'
'matched_HCs/pca12_meg14_0235_cc420486/fmraedfffpca12_meg14_0235_cc420486.nii'
'matched_HCs/pca13_meg13_0541_cc410390/fmraedfffpca13_meg13_0541_cc410390.nii'
'matched_HCs/pca14_meg12_0031/fmraedfffpca14_meg12_0031.nii'
'matched_HCs/pca15_meg14_0016_cc610405/fmraedfffpca15_meg14_0016_cc610405.nii'

'matched_HCs/pnfa1_meg14_0215_cc620572/fmraedfffpnfa1_meg14_0215_cc620572.nii'
'matched_HCs/pnfa2_meg14_0011_cc420162/fmraedfffpnfa2_meg14_0011_cc420162.nii'
'matched_HCs/pnfa3_meg14_0086_cc520055/fmraedfffpnfa3_meg14_0086_cc520055.nii'
'matched_HCs/pnfa4_meg13_0404_cc620592/fmraedfffpnfa4_meg13_0404_cc620592.nii'
'matched_HCs/pnfa5_meg08_0273/fmraedfffpnfa5_meg08_0273.nii'
'matched_HCs/pnfa6_meg14_0120_cc610061/fmraedfffpnfa6_meg14_0120_cc610061.nii'
'matched_HCs/pnfa7_meg13_0418_cc610288/fmraedfffpnfa7_meg13_0418_cc610288.nii'
'matched_HCs/pnfa8_meg14_0428_cc721291/fmraedfffpnfa8_meg14_0428_cc721291.nii'
'matched_HCs/pnfa9_meg14_0441_cc610932/fmraedfffpnfa9_meg14_0441_cc610932.nii'
'matched_HCs/pnfa10_meg08_0378/fmraedfffpnfa10_meg08_0378.nii'
};
%% MRIs
mri_path = '/imaging/rowe/archive/users/hp02/pnfa_mmn/mri_scans/';

mri_matched_hcs = {
  
'matched_HCs/bvftd1_meg08_0379/sbvftd1_meg08_0379.nii';
'matched_HCs/bvftd2_meg14_0046_cc520127/sbvftd2_meg14_0046_cc520127.nii'
'matched_HCs/bvftd3_meg13_0363_cc510629/sbvftd3_meg13_0363_cc510629.nii'
'matched_HCs/bvftd4_meg13_0479_cc420383/sbvftd4_meg13_0479_cc420383.nii'
'matched_HCs/bvftd5_meg13_0411_cc410084/sbvftd5_meg13_0411_cc410084.nii'
'matched_HCs/bvftd6_meg14_0094/sbvftd6_meg14_0094.nii'
'matched_HCs/bvftd7_meg13_0411_cc510648/sbvftd7_meg13_0411_cc510648.nii'
'matched_HCs/bvftd8_meg14_0051_cc420100/sbvftd8_meg14_0051_cc420100.nii'
'matched_HCs/bvftd9_meg14_0067_cc420157/sbvftd9_meg14_0067_cc420157.nii'
'matched_HCs/bvftd10_meg14_0114_vc4/sbvftd10_meg14_0114_vc4.nii'
'matched_HCs/bvftd11_meg13_0506_cc410179/sbvftd11_meg13_0506_cc410179.nii'
'matched_HCs/bvftd12_meg13_0403_cc510639/sbvftd12_meg13_0403_cc510639.nii'
'matched_HCs/bvftd13_meg13_0393_cc520065/sbvftd13_meg13_0393_cc520065.nii'
'matched_HCs/bvftd14_meg13_0480_cc420261/sbvftd14_meg13_0480_cc420261.nii'
'matched_HCs/bvftd15_meg13_0338_cc412021/sbvftd15_meg13_0338_cc412021.nii'
'matched_HCs/bvftd16_meg08_0368/sbvftd16_meg08_0368.nii'
'matched_HCs/bvftd17_meg12_0040/sbvftd17_meg12_0040.nii'
'matched_HCs/bvftd18_meg13_0317_cc520745/sbvftd18_meg13_0317_cc520745.nii'
'matched_HCs/bvftd19_meg14_0060_cc520253/sbvftd19_meg14_0060_cc520253.nii'
'matched_HCs/bvftd20_meg13_0401_cc420198/sbvftd20_meg13_0401_cc420198.nii'
'matched_HCs/bvftd21_meg14_0096_vc3/sbvftd21_meg14_0096_vc3.nii'
'matched_HCs/bvftd22_meg14_0271_cc520097/sbvftd22_meg14_0271_cc520097.nii'
'matched_HCs/bvftd23_meg14_0171_cc610496/sbvftd23_meg14_0171_cc610496.nii'

'matched_HCs/pca1_meg14_0044_cc510355/spca1_meg14_0044_cc510355.nii'
'matched_HCs/pca2_meg09_0087/spca2_meg09_0087.nii'
'matched_HCs/pca3_meg08_0365/spca3_meg08_0365.nii'
'matched_HCs/pca4_meg12_0520/spca4_meg12_0520.nii'
'matched_HCs/pca5_meg12_0036/spca5_meg12_0036.nii'
'matched_HCs/pca6_meg14_0056_cc410284/spca6_meg14_0056_cc410284/.nii'
'matched_HCs/pca7_meg14_0009_cc420435/spca7_meg14_0009_cc420435.nii'
'matched_HCs/pca8_meg14_0249_cc510259/spca8_meg14_0249_cc510259.nii'
'matched_HCs/pca9_meg14_0048_cc512003/spca9_meg14_0048_cc512003.nii'
'matched_HCs/pca10_meg11_0267/spca10_meg11_0267.nii'
'matched_HCs/pca11_meg13_0372_cc710446/spca11_meg13_0372_cc710446.nii'
'matched_HCs/pca12_meg14_0235_cc420486/spca12_meg14_0235_cc420486.nii'
'matched_HCs/pca13_meg13_0541_cc410390/spca13_meg13_0541_cc410390.nii'
'matched_HCs/pca14_meg12_0031/spca14_meg12_0031.nii'
'matched_HCs/pca15_meg14_0016_cc610405/spca15_meg14_0016_cc610405.nii'

'matched_HCs/pnfa1_meg14_0215_cc620572/spnfa1_meg14_0215_cc620572.nii'
'matched_HCs/pnfa2_meg14_0011_cc420162/spnfa2_meg14_0011_cc420162.nii'
'matched_HCs/pnfa3_meg14_0086_cc520055/spnfa3_meg14_0086_cc520055.nii'
'matched_HCs/pnfa4_meg13_0404_cc620592/spnfa4_meg13_0404_cc620592.nii'
'matched_HCs/pnfa5_meg08_0273/spnfa5_meg08_0273.nii'
'matched_HCs/pnfa6_meg14_0120_cc610061/spnfa6_meg14_0120_cc610061.nii'
'matched_HCs/pnfa7_meg13_0418_cc610288/spnfa7_meg13_0418_cc610288.nii'
'matched_HCs/pnfa8_meg14_0428_cc721291/spnfa8_meg14_0428_cc721291.nii'
'matched_HCs/pnfa9_meg14_0441_cc610932/spnfa9_meg14_0441_cc610932.nii'
'matched_HCs/pnfa10_meg08_0378/spnfa10_meg08_0378.nii'
};
% % % bvFTD
% % '/imaging/lh01/FTD2010/MRIscans/CBU101120/structurals/sCBU101120-0002-00001-000192-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU140288/structurals/sMR13010_CC520127-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU130867/structurals/sMR13010_CC510629-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU140097/structurals/sMR13010_CC420383-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU130971/structurals/sMR13010_CC410084-0003-00001-000192-01.nii';
% % '/imaging/lh01/FTD2010/MRIscans/CBU121067/structurals/sCBU121067-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU140091/structurals/sMR13010_CC510648-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU140144/structurals/sMR13010_CC420100-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU140242/structurals/sMR13010_CC420157-0003-00001-000192-01.nii';
% % '/imaging/rowe/archive/users/hp02/pnfa_mmn/mri_scans/vespa/meg14_0114_vc4/vc4_Structural.nii';
% % '/imaging/camcan/cc700/mri/pipeline/release003/data/aamod_copystructural_00001/CBU111038/structurals/sMR10033_CC410179-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU130803/structurals/sMR13010_CC510639-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU140053/structurals/sMR13010_CC520065-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU130918/structurals/sMR13010_CC420261-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU130811/structurals/sMR13010_CC412021-0003-00001-000192-01.nii';
% % '/imaging/lh01/FTD2010/MRIscans/CBU070612/structurals/sCBU070612-0002-00001-000160-01.nii';
% % '/imaging/lh01/FTD2010/MRIscans/CBU110542/structurals/sCBU110542-0002-00001-000192-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU130889/structurals/sMR13010_CC520745-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU140223/structurals/sMR13010_CC520253-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU131013/structurals/sMR13010_CC420198-0003-00001-000192-01.nii';
% % '/imaging/rowe/archive/users/hp02/pnfa_mmn/mri_scans/vespa/meg14_0096_vc3/vc3_Structural.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU140663/structurals/sMR13010_CC520097-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU140259/structurals/sMR13010_CC610496-0003-00001-000192-01.nii';
% % 
% % %PCA
% % '/imaging/camcan/cc700/mri/pipeline/release003/data/aamod_copystructural_00001/CBU120033/structurals/sMR10033_CC510355-0003-00001-000192-01.nii';
% % '/imaging/lh01/FTD2010/MRIscans/CBU080331/structurals/sCBU080331-0002-00001-000160-01.nii';
% % '/imaging/lh01/FTD2010/MRIscans/CBU080728/structurals/sCBU080728-0002-00001-000160-01.nii';
% % '/imaging/lh01/FTD2010/MRIscans/CBU121103/structurals/GNG/sCBU121103-0002-00001-000192-01.nii';
% % '/imaging/lh01/FTD2010/MRIscans/H17/GNG/H17_19861_mprage.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU140572/structurals/sMR13010_CC410284-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc700/mri/pipeline/release003/data/aamod_copystructural_00001/CBU121133/structurals/sMR10033_CC420435-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU140547/structurals/sMR13010_CC510259-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU140666/structurals/sMR13010_CC512003-0003-00001-000192-01.nii';
% % '/imaging/lh01/FTD2010/MRIscans/CBU111093/structurals/sCBU111093-0002-00001-000192-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU130875/structurals/sMR13010_CC710446-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU140570/structurals/sMR13010_CC420286-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc700/mri/pipeline/release003/data/aamod_copystructural_00001/CBU120290/structurals/sMR10033_CC410390-0003-00001-000192-01.nii';
% % '/imaging/lh01/FTD2010/MRIscans/CBU120153/structurals/sCBU120153-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU140149/structurals/sMR13010_CC610405-0003-00001-000097-01.nii';
% % 
% % % PNFA
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU140612/structurals/sMR13010_CC620572-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU140132/structurals/sMR13010_CC420162-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU140279/structurals/sMR13010_CC520055-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU130933/structurals/sMR13010_CC620592-0003-00001-000192-01.nii';
% % '/imaging/lh01/FTD2010/MRIscans/CBU070619/structurals/y_sCBU070619-0002-00001-000160-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU140124/structurals/sMR13010_CC610061-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc280/mri/pipeline/release003/Structural/aamod_copystructural_00001/CBU130932/structurals/sMR13010_CC610288-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc700/mri/pipeline/release003/data/aamod_copystructural_00001/CBU130005/structurals/sMR10033_CC721291-0003-00001-000192-01.nii';
% % '/imaging/camcan/cc700/mri/pipeline/release003/data/aamod_copystructural_00001/CBU111154/structurals/sMR10033_CC610392-0003-00001-000192-01.nii';
% % '/imaging/lh01/FTD2010/MRIscans/CBU080339/structurals/y_sCBU080339-0002-00001-000160-01.nii';
% % };
