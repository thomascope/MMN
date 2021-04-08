function extDCM_nodetest(dirname_DCM,filestem,conditions,groups,p,all_names)
% A function to be called in parallel to do firstlevel DCMs, for example,
% with the existing structures p and all_names from ICA_denoise_Hollydata mainscript
% dirname_DCM = '/imaging/mlr/users/tc02/Holly_MMN/extDCMs/';
% filestem = 'b8LFP_s_-100_500_LOR_fmbraedfffM';
% conditions = {'STD', 'DVT'};
% all_combinations = combvec(groups,1:length(conditions));
% parfor this_comb = 1:length(all_combinations)
%     k = all_combinations(1,this_comb)
%     c = all_combinations(2,this_comb)
%     extDCM_firstlevel_PEB(dirname_DCM,conditions(c),k,p,all_names)
% end

%% PEB 1 - THE FULL MODEL:

% this is repeated for each set of DCMs you did. For instance, I had the
% controls and patients for the standard and deviant conditions in
% different folders, making 4 repitions of this section. The full model is
% reduced for the different folders separately, meaning there is no
% group-mean assumption across the conditions - only across the
% participants in each group.

% First, go to the directory where you saved one of your DCMs and create
% your file list to work with:


for k = groups
    for c = 1:length(conditions)
                
        this_group = find(p.group==k);
        dcm_files={};
        for subj = 1:length(this_group)
            if exist([dirname_DCM filestem all_names{this_group(subj)} '_dcm_' conditions{c} '.mat'],'file')
                dcm_files{end+1} = [dirname_DCM filestem all_names{this_group(subj)} '_dcm_' conditions{c} '.mat'];
            else
                warning([dirname_DCM filestem all_names{this_group(subj)} '_dcm_' conditions{c} '.mat does not exist'])
            end
        end
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
        % output, not the PEB. Note that this takes the place of the more usual spm_dcm_peb_fit as we do not want to re-estimate the DCM:
        % cd to a dummy directory as, if parallelised, tmp.mat files can
        % overwrite each other
        mkdir([pwd filesep 'nodetest_folders' filesep 'tempdir_' p.diagnosis_list{k} '_' conditions{c}])
        cd([pwd filesep 'nodetest_folders' filesep 'tempdir_' p.diagnosis_list{k} '_' conditions{c}])
        [PEB,DCM] = spm_dcm_peb(DCM,M,'all');
        DCM = spm_dcm_reduce(DCM,rE,rC);
        cd(this_dir)
    end
end



end