function extDCM_overall_PEB(dirname_DCM,filestem,conditions,groups,p,all_names,k)

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

if k == 1
    PEB_ABC = load('PEB_ABC.mat','PEB_ABC')
    PEB_ABC = PEB_ABC.PEB_ABC;
    X = ones(length(groups),1); %Group mean
    for i = 1:length(groups) %Every pairwise contrast once
        j = i;
        while j<=length(groups)-1
            j = j+1;
            X(i,end+1) = 1;
            X(j,end) = -1;
        end
    end
    
    mkdir([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_ABC_' char(join(conditions,'_')) '_Overall'])
    cd([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_ABC_' char(join(conditions,'_')) '_Overall'])
    [PEB_ABC_Overall,DCM_ABC_Overall] = spm_dcm_peb(PEB_ABC,X,{'A','B','C'});
    cd([dirname_DCM 'PEB_secondlevel'])
    save('PEB_ABC_Overall.mat','PEB_ABC_Overall','DCM_ABC_Overall')
elseif k == 2
    PEB_HD = load('PEB_HD.mat','PEB_HD')
    PEB_HD = PEB_HD.PEB_HD;
    X = ones(length(groups),1); %Group mean
    for i = 1:length(groups) %Every pairwise contrast once
        j = i;
        while j<=length(groups)-1
            j = j+1;
            X(i,end+1) = 1;
            X(j,end) = -1;
        end
    end
    
    mkdir([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_HD_' char(join(conditions,'_')) '_Overall'])
    cd([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_HD_' char(join(conditions,'_')) '_Overall'])
    [PEB_HD_Overall,DCM_HD_Overall] = spm_dcm_peb(PEB_HD,X,{'H','D'});
    cd([dirname_DCM 'PEB_secondlevel'])
    save('PEB_HD_Overall.mat','PEB_HD_Overall','DCM_HD_Overall')
end


if k == 1
    cd([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_ABC_' char(join(conditions,'_')) '_Overall'])
    % find the posterior p-vals for the parameters in this PEB
    BMA_ABC_Overall = spm_dcm_peb_bmc(PEB_ABC_Overall);
    BMA_ABC_Overall.Ep = reshape(BMA_ABC_Overall.Ep,[],size(X,2));
    BMA_ABC_Overall.Cp = reshape(BMA_ABC_Overall.Cp,[],size(X,2));
    BMA_ABC_Overall.Pp = reshape(BMA_ABC_Overall.Pp,[],size(X,2));
    cd([dirname_DCM 'PEB_secondlevel'])
    save('PEB_ABC_Overall.mat','BMA_ABC_Overall','-append')
elseif k == 2
    cd([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_HD_' char(join(conditions,'_')) '_Overall'])
    BMA_HD_Overall = spm_dcm_peb_bmc(PEB_HD_Overall);
    BMA_HD_Overall.Ep = reshape(BMA_HD_Overall.Ep,[],size(X,2));
    BMA_HD_Overall.Cp = reshape(BMA_HD_Overall.Cp,[],size(X,2));
    BMA_HD_Overall.Pp = reshape(BMA_HD_Overall.Pp,[],size(X,2));
    cd([dirname_DCM 'PEB_secondlevel'])
    save('PEB_HD_Overall.mat','BMA_HD_Overall','-append')
end

    function restore_env(old_path)
        path(old_path);
    end

end

