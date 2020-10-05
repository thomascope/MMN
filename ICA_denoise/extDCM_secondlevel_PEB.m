function extDCM_secondlevel_PEB(dirname_DCM,filestem,conditions,groups,p,all_names)

%% PEB 2 - the second level design matrix

% load all the reduced DCMs created in the above step into a single
% structure (option here of doing a single shot or a group-wise PEB of PEBs, see https://en.wikibooks.org/wiki/SPM/Parametric_Empirical_Bayes_(PEB)#Estimate_a_second_level_PEB_(Parametric_Empirical_Bayes)_model):
% GCM = {};
% X1 = [];
% i = 0;
% for k = groups
%     for c = 1:length(conditions)
%         i = i+1;
%         disp(['loading file ' num2str(i) ' of ' num2str(length(groups)*length(conditions))])
%         GCMtmp = load([dirname_DCM 'PEB_firstlevel' filesep 'First_level_PEB_' p.diagnosis_list{k} '_' conditions{c} '.mat']);
%         GCMtmp = GCMtmp.DCM;
%         GCM(end+1:end+length(GCMtmp),1) = GCMtmp;
%         X1(end+1:end+length(GCMtmp),end+1) = ones(length(GCMtmp),1);
%     end
% end
% 
GCM = {};
X1 = {};
i = 0;
for k = groups
    for c = 1:length(conditions)
        i = i+1;
        disp(['loading file ' num2str(i) ' of ' num2str(length(groups)*length(conditions))])
        GCMtmp = load([dirname_DCM 'PEB_firstlevel' filesep 'First_level_PEB_' p.diagnosis_list{k} '_' conditions{c} '.mat']);
        GCMtmp = GCMtmp.DCM;
        try
            GCM{k}(end+1:end+length(GCMtmp),1) = GCMtmp;
            X1{k}(end+1:end+length(GCMtmp),end+1) = ones(length(GCMtmp),1);
        catch
            GCM{k} = {};
            X1{k} = [];
            GCM{k}(end+1:end+length(GCMtmp),1) = GCMtmp;
            X1{k}(end+1:end+length(GCMtmp),end+1) = ones(length(GCMtmp),1);
        end
    end
end

if iscell(X1)
    PEB_ABC2 = cell(2*length(groups),1);
    PEB_HD2 = cell(2*length(groups),1);
    DCM_ABC2 = cell(2*length(groups),1);
    DCM_HD2 = cell(2*length(groups),1);
    BMA_ABC2 = cell(2*length(groups),1);
    BMA_HD2 = cell(2*length(groups),1);

    parfor k = 1:2*length(groups)
        
        % Second level, per group, PEB
        % Make a model structure, M, with a design matrix specified in .X. The
        % matrix needs a column of ones first, which represents the group mean. In
        % this example the rest of the columns represent the condition differences,
        % (from the differen DCM_all folders) with an example interaction in the
        % final column:
        
        %warning on verbose
        warning('error','MATLAB:nearlySingularMatrix') % Make the PEB fail if enters a singular matrix state
        
        if k <=length(groups)
            X = ones(size(X1{k},1),1);
            for i = 1:size(p.contrast_weights,1)
                X(:,end+1) = sum(p.contrast_weights(i,:).*X1{k},2);
            end
            % figure;imagesc(X);colorbar %Visualise design matrix for sanity check
            M = GCM{k}{1}.M;
            M.X = X;
            mkdir([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_ABC2' p.diagnosis_list{k}])
            cd([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_ABC2' p.diagnosis_list{k}])
            
            try
                [PEB_ABC2{k},DCM_ABC2{k}] = spm_dcm_peb(GCM{k},M,{'A','B','C'}) % Intrinsic and Extrinsic connections, plus input, plus inter-regional delays, overall
            catch
%                 try
%                     [PEB_ABC2{k},DCM_ABC2{k}] = spm_dcm_peb(GCM{k},M,{'A','B','C','D'}) % Add an extra parameter of no interest if there are optimisation errors
%                 catch
%                     error(['The failed condition was ABC2 group ' num2str(k)])
%                 end
            end
            % find the posterior p-vals for the parameters in this PEB
            BMA_ABC2{k} = spm_dcm_peb_bmc(PEB_ABC2{k});
            BMA_ABC2{k}.Ep = reshape(BMA_ABC2{k}.Ep,[],size(M.X,2));
            BMA_ABC2{k}.Cp = reshape(BMA_ABC2{k}.Cp,[],size(M.X,2));
            BMA_ABC2{k}.Pp = reshape(BMA_ABC2{k}.Pp,[],size(M.X,2));
        else
            X = ones(size(X1{k-length(groups)},1),1);
            for i = 1:size(p.contrast_weights,1)
                X(:,end+1) = sum(p.contrast_weights(i,:).*X1{k-length(groups)},2);
            end
            % figure;imagesc(X);colorbar %Visualise design matrix for sanity check
            M = GCM{k-length(groups)}{1}.M;
            M.X = X;
            mkdir([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_HD2' p.diagnosis_list{k-length(groups)}])
            cd([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_HD2' p.diagnosis_list{k-length(groups)}])
            
            try
                [PEB_HD2{k},DCM_HD2{k}] = spm_dcm_peb(GCM{k-length(groups)},M,{'H','D'}) % Intrinsic connections by population
            catch
%                 try
%                     [PEB_HD2{k},DCM_HD2{k}] = spm_dcm_peb(GCM{k-length(groups)},M,{'H','D'}) % Add an extra parameter of no interest if there are optimisation errors
%                 catch
%                     error(['The failed condition was H group ' num2str(k-length(groups)) ', k = ' num2str(k)])
%                 end
            end
            % find the posterior p-vals for the parameters in this PEB
            BMA_HD2{k} = spm_dcm_peb_bmc(PEB_HD2{k});
            BMA_HD2{k}.Ep = reshape(BMA_HD2{k}.Ep,[],size(M.X,2));
            BMA_HD2{k}.Cp = reshape(BMA_HD2{k}.Cp,[],size(M.X,2));
            BMA_HD2{k}.Pp = reshape(BMA_HD2{k}.Pp,[],size(M.X,2));
        end
    end
    for k = 1:length(groups)
        PEB_HD2(1) = [];
        DCM_HD2(1) = [];
        BMA_HD2(1) = [];
        PEB_ABC2(end) = [];
        DCM_ABC2(end) = [];
        BMA_ABC2(end) = [];
    end
    cd([dirname_DCM 'PEB_secondlevel'])
    save('PEB_ABC2.mat','PEB_ABC2','DCM_ABC2','BMA_ABC2')
    save('PEB_HD2.mat','PEB_HD2','DCM_HD2','BMA_HD2')
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
            mkdir([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_ABC2_Overall'])
            cd([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_ABC2_Overall'])
            [PEB_ABC2_Overall,DCM_ABC2_Overall] = spm_dcm_peb(PEB_ABC2,X)
            % find the posterior p-vals for the parameters in this PEB
            BMA_ABC2_Overall = spm_dcm_peb_bmc(PEB_ABC2_Overall);
            BMA_ABC2_Overall.Ep = reshape(BMA_ABC2_Overall.Ep,[],size(M.X,2));
            BMA_ABC2_Overall.Cp = reshape(BMA_ABC2_Overall.Cp,[],size(M.X,2));
            BMA_ABC2_Overall.Pp = reshape(BMA_ABC2_Overall.Pp,[],size(M.X,2));
        else
            mkdir([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_HD2_Overall'])
            cd([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_HD2_Overall'])
            [PEB_HD2_Overall,DCM_HD2_Overall] = spm_dcm_peb(PEB_HD2,X)
            % find the posterior p-vals for the parameters in this PEB
            BMA_HD2_Overall = spm_dcm_peb_bmc(PEB_HD2_Overall);
            BMA_HD2_Overall.Ep = reshape(BMA_HD2_Overall.Ep,[],size(M.X,2));
            BMA_HD2_Overall.Cp = reshape(BMA_HD2_Overall.Cp,[],size(M.X,2));
            BMA_HD2_Overall.Pp = reshape(BMA_HD2_Overall.Pp,[],size(M.X,2));
        end
    end
    cd([dirname_DCM 'PEB_secondlevel'])
    save('PEB_ABC2_Overall.mat','PEB_ABC2_Overall','DCM_ABC2_Overall','BMA_ABC2_Overall')
    save('PEB_HD2_Overall.mat','PEB_HD2_Overall','DCM_HD2_Overall','BMA_HD2_Overall')
    delete(gcp)
else
    error('I have only written a two-stage implementation, feel free to add a one shot version here!')
end
% 
% % Now you have your DCMs and a design matrix, you need to decide which
% % parameters you a interested in. This is basically a list of connections, 
% % either extrinsic or intrinsic, that you are interested in. e.g. It could 
% % just be all the extrinsics: ... 
% params = {'A'};
% % ... or it could be a specific list of intrinsics - in the below example 
% % all GABA-A synapses:
% params = {'H(1,1,1,3)'
%           'H(2,2,1,3)'
%           'H(3,3,1,3)'
%           'H(4,4,1,3)'
%           'H(5,5,1,3)'
%           'H(6,6,1,3)'
%           'H(2,3,1,3)'
%           'H(1,3,1,3)'
%           'H(4,5,1,3)'
%           'H(6,5,1,3)'};
% 
% % Run the second level PEB and save the results. The .Ep field in the PEB
% % has a column for each column of the design matrix:
% [PEB_2,DCM_2] = spm_dcm_peb(GCM,M,params);
% save('PEB_2','PEB_2','DCM_2')
% 
% % The above results can then be used to find the posterior p-vals for the
% % parameters in this PEB (found in .pP). The relevant fields need reshaping
% % to the size of the design matrix:
% BMA_2 = spm_dcm_peb_bmc(PEB_2);
% BMA_2.Ep = reshape(BMA_2.Ep,[],size(M.X,2));
% BMA_2.Cp = reshape(BMA_2.Cp,[],size(M.X,2));
% BMA_2.Pp = reshape(BMA_2.Pp,[],size(M.X,2));
% save('PEB_2','BMA_2','-append')
% 
% % You are now ready to view some results. I have created a plotting script 
% % based on my 6-node, 6-population model. It shows the design matrix and 
% % then a row of graph plots, one plot for each column of the design matrix 
% % (not including the first column, as this is how the PEB has adjusted the 
% % parameters for the group mean - not the constrasts you are interested in).
% % To use it for a single BMA to view both the parameter estimates and their
% % associated P-value:
% clear B
% B.BMA.bma = BMA;
% plotPEBBMA(B,'Pp','cax',.5)
% 
% % To use it for multiple BMAs in one figure (only allowing you to look at 
% % the parameter estimates, not their P-values):
% clear B
% B.BMA1.bma = BMA1;
% B.BMA2.bma = BMA2;
% plotPEBBMA(B,'Ep','cax',.5)
% 
