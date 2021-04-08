function node_test_script

spmpath = '/group/language/data/thomascope/MMN/ICA_denoise/Tallie_extDCM/spm12_latestTA/';
thisspm = which('spm');
if ~strcmp(thisspm(1:end-5), spmpath)
    rmpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/'));
    rmpath(genpath('/group/language/data/thomascope/spm12_fil_r6906/'));
    addpath(spmpath)
    spm eeg
end

addpath('/group/language/data/thomascope/MMN/ICA_denoise/Tallie_extDCM/mfiles_also_needed');

PEB_HD = load('/imaging/mlr/users/tc02/Holly_MMN/extDCMs_hE6/PEB_secondlevel/PEB_HD.mat','PEB_HD');
PEB_HD = PEB_HD.PEB_HD;

groups = 1:5;
X = ones(length(groups),1); %Group mean
for i = 1:length(groups) %Every pairwise contrast once
    j = i;
    while j<=length(groups)-1
        j = j+1;
        X(i,end+1) = 1;
        X(j,end) = -1;
    end
end

mkdir([pwd '/node_test_tempdir'])
cd([pwd '/node_test_tempdir'])
[PEB_HD_Overall,DCM_HD_Overall] = spm_dcm_peb(PEB_HD,X,{'H','D'})
% find the posterior p-vals for the parameters in this PEB
BMA_HD_Overall = spm_dcm_peb_bmc(PEB_HD_Overall);
