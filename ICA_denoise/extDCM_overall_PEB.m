function extDCM_overall_PEB(dirname_DCM,filestem,conditions,groups,p,all_names)

old_path = path;
cleanupObj = onCleanup(@()restore_env(old_path));

spmpath = '/group/language/data/thomascope/MMN/ICA_denoise/Tallie_extDCM/spm12_latestTA/';
thisspm = which('spm');
if ~strcmp(thisspm(1:end-5), spmpath)
    rmpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/'));
    rmpath(genpath('/group/language/data/thomascope/spm12_fil_r6906/'));
    addpath(spmpath)
    spm eeg
end

addpath('/group/language/data/thomascope/MMN/ICA_denoise/Tallie_extDCM/mfiles_also_needed');

cd([dirname_DCM 'PEB_secondlevel'])
load('temp_workspace.mat')
load('PEB_ABC.mat','PEB_ABC','DCM_ABC','BMA_ABC')
load('PEB_HD.mat','PEB_HD','DCM_HD','BMA_HD')

parfor k = 1:2
    if k == 1
        M = PEB_ABC{1}.M
        X = ones(length(groups),1); %Group mean
        for i = 1:length(groups) %Every pairwise contrast once
            j = i;
            while j<=length(groups)-1
                j = j+1;
                X(i,end+1) = 1;
                X(j,end) = -1;
            end
        end
        M.X = X;
        
        mkdir([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_ABC_' char(join(conditions,'_')) '_Overall'])
        cd([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_ABC_' char(join(conditions,'_')) '_Overall'])
        [PEB_ABC_Overall,DCM_ABC_Overall] = spm_dcm_peb(PEB_ABC,M,{'A','B','C'})
        % find the posterior p-vals for the parameters in this PEB
        BMA_ABC_Overall = spm_dcm_peb_bmc(PEB_ABC_Overall);
        BMA_ABC_Overall.Ep = reshape(BMA_ABC_Overall.Ep,[],size(M.X,2));
        BMA_ABC_Overall.Cp = reshape(BMA_ABC_Overall.Cp,[],size(M.X,2));
        BMA_ABC_Overall.Pp = reshape(BMA_ABC_Overall.Pp,[],size(M.X,2));
    else
        M = PEB_HD{1}.M
        X = ones(length(groups),1); %Group mean
        for i = 1:length(groups) %Every pairwise contrast once
            j = i;
            while j<=length(groups)-1
                j = j+1;
                X(i,end+1) = 1;
                X(j,end) = -1;
            end
        end
        M.X = X;
        
        mkdir([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_HD_' char(join(conditions,'_')) '_Overall'])
        cd([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_HD_' char(join(conditions,'_')) '_Overall'])
        [PEB_HD_Overall,DCM_HD_Overall] = spm_dcm_peb(PEB_HD,M,{'H','D'})
        % find the posterior p-vals for the parameters in this PEB
        BMA_HD_Overall = spm_dcm_peb_bmc(PEB_HD_Overall);
        BMA_HD_Overall.Ep = reshape(BMA_HD_Overall.Ep,[],size(M.X,2));
        BMA_HD_Overall.Cp = reshape(BMA_HD_Overall.Cp,[],size(M.X,2));
        BMA_HD_Overall.Pp = reshape(BMA_HD_Overall.Pp,[],size(M.X,2));
    end
end
cd([dirname_DCM 'PEB_secondlevel'])
save('PEB_ABC_Overall.mat','PEB_ABC_Overall','DCM_ABC_Overall','BMA_ABC_Overall')
save('PEB_HD_Overall.mat','PEB_HD_Overall','DCM_HD_Overall','BMA_HD_Overall')

    function restore_env(old_path)
        path(old_path);
    end

end

