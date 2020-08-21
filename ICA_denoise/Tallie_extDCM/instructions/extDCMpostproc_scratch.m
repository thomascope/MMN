
%
% extDCMpostproc_notes
%


%% PEB 1 - THE FULL MODEL:

% this is repeated for each set of DCMs you did. For instance, I had the
% controls and patients for the standard and deviant conditions in
% different folders, making 4 repitions of this section. The full model is
% reduced for the different folders separately, meaning there is no
% group-mean assumption across the conditions - only across the
% participants in each group.

% First, go to the directory where you saved one of your DCMs and create
% your file list to work with:

dirname_DCM = '/imaging/tc02/Holly_MMN/extDCMs/';
conditions = {'STD', 'DVT'};
all_combinations = combvec(1:max(p.group),1:length(conditions));
parfor this_comb = 1:length(all_combinations)
    k = all_combinations(1,this_comb)
    c = all_combinations(2,this_comb)
% for k = 1:max(p.group)
%     for c = 1:length(conditions)
        this_group = find(p.group==k);
        dcm_files={};
        for subj = 1:length(this_group)
            if exist([dirname_DCM 'b8LFP_s_-100_500_LOR_fmbraedfffM' all_names{this_group(subj)} '_dcm_' conditions{c} '.mat'],'file')
            dcm_files{end+1} = [dirname_DCM 'b8LFP_s_-100_500_LOR_fmbraedfffM' all_names{this_group(subj)} '_dcm_' conditions{c} '.mat'];
            else
                warning([dirname_DCM 'b8LFP_s_-100_500_LOR_fmbraedfffM' all_names{this_group(subj)} '_dcm_' conditions{c} '.mat does not exist'])
            end
        end

    save('dcm_files','dcm_files')

    %  Load the DCMs into a cell array and set up the variables for the PEB:
    clear DCM
    for i = 1:length(dcm_files)
        DCM{i,1} = load(dcm_files{i});
        DCM{i,1} = DCM{i,1}.DCM;
    end
    [~,rC,rE] = spm_find_pC(DCM{1});
    clear M
    M.Q = 'single';
    M.bE = rE;
    M.bC = rC;

    % run the PEB for all fields, giving the group-level result in PEB and the
    % reduced values for all individuals in DCM (N.B. only .Ep is updated in
    % the new DCM structure). In this case we are interested in the DCM
    % output, not the PEB:
    [PEB,DCM] = spm_dcm_peb(DCM,M,'all');
    DCM = spm_dcm_reduce(DCM,rE,rC);
    
    if ~exist([dirname_DCM 'PEB_firstlevel' filesep])
mkdir([dirname_DCM 'PEB_firstlevel' filesep])
    end

    save([dirname_DCM 'PEB_firstlevel' filesep 'First_level_PEB_' p.diagnosis_list{k} '_' conditions{c} '.mat'],'DCM','PEB')
end
end

%% PEB 2 - the second level design matrix

% load all the reduced DCMs created in the above step into a single
% structure:
GCM = {};
X1 = [];
for k = 1:length(subdirname_DCM)
    GCMtmp = load([dirname_DCM subdirname_DCM{k} nam]);
    GCMtmp = GCMtmp.DCM;
    GCM(end+1:end+length(GCMtmp),1) = GCMtmp;
    X1(end+1:end+length(GCMtmp),end+1) = ones(length(GCMtmp),1);
end

% Make a model structure, M, with a design matrix specified in .X. The
% matrix needs a column of ones first, which represents the group mean. In
% this example the rest of the columns represent the condition differences,
% (from the differen DCM_all folders) with an example interaction in the
% final column:
X = [X1(:,1)-X1(:,2) X1(:,3)-X1(:,4)]; % condition differences
X(:,4) = X(:,2).*X(:,3); % interaction
M = GCM{1}.M;
M.X = [ones(size(GCM,1),1) X];
% N.B. the above is just an illustrative example. In my own data,
% conditions had to be further separated within these groups. Be aware you
% will need to build your own design matrix, X, for your individual needs.

% Now you have your DCMs and a design matrix, you need to decide which
% parameters you a interested in. This is basically a list of connections, 
% either extrinsic or intrinsic, that you are interested in. e.g. It could 
% just be all the extrinsics: ... 
params = {'A'};
% ... or it could be a specific list of intrinsics - in the below example 
% all GABA-A synapses:
params = {'H(1,1,1,3)'
          'H(2,2,1,3)'
          'H(3,3,1,3)'
          'H(4,4,1,3)'
          'H(5,5,1,3)'
          'H(6,6,1,3)'
          'H(2,3,1,3)'
          'H(1,3,1,3)'
          'H(4,5,1,3)'
          'H(6,5,1,3)'};

% Run the second level PEB and save the results. The .Ep field in the PEB
% has a column for each column of the design matrix:
[PEB_2,DCM_2] = spm_dcm_peb(GCM,M,params);
save('PEB_2','PEB_2','DCM_2')

% The above results can then be used to find the posterior p-vals for the
% parameters in this PEB (found in .pP). The relevant fields need reshaping
% to the size of the design matrix:
BMA_2 = spm_dcm_peb_bmc(PEB_2);
BMA_2.Ep = reshape(BMA_2.Ep,[],size(M.X,2));
BMA_2.Cp = reshape(BMA_2.Cp,[],size(M.X,2));
BMA_2.Pp = reshape(BMA_2.Pp,[],size(M.X,2));
save('PEB_2','BMA_2','-append')

% You are now ready to view some results. I have created a plotting script 
% based on my 6-node, 6-population model. It shows the design matrix and 
% then a row of graph plots, one plot for each column of the design matrix 
% (not including the first column, as this is how the PEB has adjusted the 
% parameters for the group mean - not the constrasts you are interested in).
% To use it for a single BMA to view both the parameter estimates and their
% associated P-value:
clear B
B.BMA.bma = BMA;
plotPEBBMA(B,'Pp','cax',.5)

% To use it for multiple BMAs in one figure (only allowing you to look at 
% the parameter estimates, not their P-values):
clear B
B.BMA1.bma = BMA1;
B.BMA2.bma = BMA2;
plotPEBBMA(B,'Ep','cax',.5)


