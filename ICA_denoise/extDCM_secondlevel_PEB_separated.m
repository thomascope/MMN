function extDCM_secondlevel_PEB_separated(dirname_DCM,conditions,groups,p,PEB_focus,regions,conductances)

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

%% PEB 2 - the second level design matrix

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

these_conditions = find(contains(p.conditions,conditions));
all_contrast_cols = 1:size(p.contrast_weights,2);
other_cols = setdiff(all_contrast_cols,these_conditions);
these_contrasts = [];
for i = 1:size(p.contrast_weights,1)
    if sum(p.contrast_weights(i,these_conditions)~=0)>=2&&all(p.contrast_weights(i,other_cols)==0)
        these_contrasts(end+1) = i;
    end
end
% Pare down the contrasts to only those in the included conditions
p.contrast_weights = p.contrast_weights(these_contrasts,these_conditions);
p.contrast_labels = p.contrast_labels(these_contrasts);

%get rid of duplicates and inverse contrasts
duplicates = [];
for i = 1:size(p.contrast_weights,1)
    for j = (i+1):size(p.contrast_weights,1)
        if all(p.contrast_weights(i,:)==-p.contrast_weights(j,:))||all(p.contrast_weights(i,:)==p.contrast_weights(j,:))
            duplicates(end+1) = j;
        end
    end
end
p.contrast_weights(duplicates,:) = [];
p.contrast_labels(duplicates) = [];

PEB = cell(length(PEB_focus),length(groups));
DCM = cell(length(PEB_focus),length(groups));
BMA = cell(length(PEB_focus),length(groups));
for i = 1:length(PEB_focus)
    for k = groups
        these_labels = regexp(PEB_focus(i),'\d*','Match');
        if isempty(these_labels{1})
            this_tail = PEB_focus{i};
        elseif length(these_labels{1})==2
            this_tail = [PEB_focus{i}(1) '_' regions{str2num(these_labels{1}{1})} '_' conductances{str2num(these_labels{1}{2})}];
        else
            error('Neither zero nor two numbers in the PEB focus - I do not know how to interpret this')
        end
        
        X = ones(size(X1{k},1),1);
        for i = 1:size(p.contrast_weights,1)
            X(:,end+1) = sum(p.contrast_weights(i,:).*X1{k},2);
        end
        % figure;imagesc(X);colorbar %Visualise design matrix for sanity check
        M = GCM{k}{1}.M;
        M.X = X;
        
        mkdir([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_' this_tail '_' char(join(conditions,'_')) '_' p.diagnosis_list{k}])
        cd([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_' this_tail '_' char(join(conditions,'_')) '_' p.diagnosis_list{k}])
        
        try
            [PEB{i,k},DCM{i,k}] = spm_dcm_peb(GCM{k},M,PEB_focus(i)); % Second level PEB
        catch
            warning(['The failed condition was ' this_tail ' group ' num2str(k)])
        end
        
        BMA{i,k} = spm_dcm_peb_bmc(PEB{i,k});
        BMA{i,k}.Ep = reshape(BMA{i,k}.Ep,[],size(M.X,2));
        BMA{i,k}.Cp = reshape(BMA{i,k}.Cp,[],size(M.X,2));
        BMA{i,k}.Pp = reshape(BMA{i,k}.Pp,[],size(M.X,2));
        
        PEB = PEB{i,k};
        DCM = DCM{i,k};
        BMA = BMA{i,k};
        cd([dirname_DCM 'PEB_secondlevel'])
        save(['PEB_' this_tail '_' p.diagnosis_list{k} '.mat'],'PEB','DCM','BMA')
        
    end
    
end

function restore_env(old_path)
path(old_path);
end

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

