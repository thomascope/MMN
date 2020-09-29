%parfor k = 1:2*length(groups)
parfor k = 1:10
    
    % Second level, per group, PEB
    % Make a model structure, M, with a design matrix specified in .X. The
    % matrix needs a column of ones first, which represents the group mean. In
    % this example the rest of the columns represent the condition differences,
    % (from the differen DCM_all folders) with an example interaction in the
    % final column:
    
    
    if k <=length(groups)
        X = ones(size(X1{k},1),1);
        for i = 1:size(p.contrast_weights,1)
            X(:,end+1) = sum(p.contrast_weights(i,:).*X1{k},2);
        end
        % figure;imagesc(X);colorbar %Visualise design matrix for sanity check
        M = GCM{k}{1}.M;
        M.X = X;
        mkdir([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_ABC' p.diagnosis_list{k}])
        cd([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_ABC' p.diagnosis_list{k}])
        
        [PEB_ABC{k},DCM_ABC{k}] = spm_dcm_peb(GCM{k},M,{'A','B','C'}) % Intrinsic and Extrinsic connections, plus input, overall
        % find the posterior p-vals for the parameters in this PEB
        BMA_ABC{k} = spm_dcm_peb_bmc(PEB_ABC{k});
        BMA_ABC{k}.Ep = reshape(BMA_ABC{k}.Ep,[],size(M.X,2));
        BMA_ABC{k}.Cp = reshape(BMA_ABC{k}.Cp,[],size(M.X,2));
        BMA_ABC{k}.Pp = reshape(BMA_ABC{k}.Pp,[],size(M.X,2));
    else
        X = ones(size(X1{k-length(groups)},1),1);
        for i = 1:size(p.contrast_weights,1)
            X(:,end+1) = sum(p.contrast_weights(i,:).*X1{k-length(groups)},2);
        end
        % figure;imagesc(X);colorbar %Visualise design matrix for sanity check
        M = GCM{k-length(groups)}{1}.M;
        M.X = X;
        mkdir([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_HD' p.diagnosis_list{k-length(groups)}])
        cd([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_HD' p.diagnosis_list{k-length(groups)}])
        [PEB_HD{k},DCM_HD{k}] = spm_dcm_peb(GCM{k-length(groups)},M,{'H','D'}) % Intrinsic connections by population and inter-regional delays
        % find the posterior p-vals for the parameters in this PEB
        BMA_HD{k} = spm_dcm_peb_bmc(PEB_HD{k});
        BMA_HD{k}.Ep = reshape(BMA_HD{k}.Ep,[],size(M.X,2));
        BMA_HD{k}.Cp = reshape(BMA_HD{k}.Cp,[],size(M.X,2));
        BMA_HD{k}.Pp = reshape(BMA_HD{k}.Pp,[],size(M.X,2));
    end
end
for k = 1:length(groups)
    PEB_HD(1) = [];
    DCM_HD(1) = [];
    BMA_HD(1) = [];
end
cd([dirname_DCM 'PEB_secondlevel'])
save('PEB_ABC.mat','PEB_ABC','DCM_ABC','BMA_ABC')
save('PEB_HD.mat','PEB_HD','DCM_HD','BMA_HD')
%Now take it to a PEB of PEBs
parfor k = 1:2
    X = ones(length(groups),1) %Group mean
    for i = 1:length(groups) %Every pairwise contrast once
        j = i
        while j<=length(groups)
            j = j+1
            X(i,end+1) = 1;
            X(j,end) = -1;
        end
    end
    if k == 1
        mkdir([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_ABC_Overall'])
        cd([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_ABC_Overall'])
        [PEB_ABC_Overall,DCM_ABC_Overall] = spm_dcm_peb(PEB_ABC,X)
        % find the posterior p-vals for the parameters in this PEB
        BMA_ABC_Overall = spm_dcm_peb_bmc(PEB_ABC_Overall);
        BMA_ABC_Overall.Ep = reshape(BMA_ABC_Overall.Ep,[],size(M.X,2));
        BMA_ABC_Overall.Cp = reshape(BMA_ABC_Overall.Cp,[],size(M.X,2));
        BMA_ABC_Overall.Pp = reshape(BMA_ABC_Overall.Pp,[],size(M.X,2));
    else
        mkdir([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_HD_Overall'])
        cd([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_HD_Overall'])
        [PEB_HD_Overall,DCM_HD_Overall] = spm_dcm_peb(PEB_HD,X)
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
delete(gcp)