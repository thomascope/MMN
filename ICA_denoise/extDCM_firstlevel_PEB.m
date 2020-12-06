function extDCM_firstlevel_PEB(dirname_DCM,filestem,conditions,groups,p,all_names)
% A function to be called in parallel to do firstlevel DCMs, for example,
% with the existing structures p and all_names from ICA_denoise_Hollydata mainscript
% dirname_DCM = '/imaging/tc02/Holly_MMN/extDCMs/';
% filestem = 'b8LFP_s_-100_500_LOR_fmbraedfffM';
% conditions = {'STD', 'DVT'};
% all_combinations = combvec(groups,1:length(conditions));
% parfor this_comb = 1:length(all_combinations)
%     k = all_combinations(1,this_comb)
%     c = all_combinations(2,this_comb)
%     extDCM_firstlevel_PEB(dirname_DCM,conditions(c),k,p,all_names)
% end

old_path = path;
cleanupObj = onCleanup(@()restore_env(old_path));
this_dir = pwd;

spmpath = '/group/language/data/thomascope/MMN/ICA_denoise/Tallie_extDCM/spm12_latestTA/';
thisspm = which('spm');
if ~strcmp(thisspm(1:end-5), spmpath)
    rmpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/'));
    rmpath(genpath('/group/language/data/thomascope/spm12_fil_r6906/'));
    addpath(spmpath)
    spm eeg
end

addpath('/group/language/data/thomascope/MMN/ICA_denoise/Tallie_extDCM/mfiles_also_needed');

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
                if exist([dirname_DCM filestem all_names{this_group(subj)} '_1_dcm_' conditions{c} '.mat'],'file')
                    try
                        movefile([dirname_DCM filestem all_names{this_group(subj)} '_1_dcm_' conditions{c} '.mat'],[dirname_DCM filestem all_names{this_group(subj)} '_dcm_' conditions{c} '.mat'])
                    catch
                        error([dirname_DCM filestem all_names{this_group(subj)} '_dcm_' conditions{c} '.mat does not exist'])
                    end
                end
            end
        end
        
        save(['dcm_files'  p.diagnosis_list{k} '_' conditions{c}],'dcm_files')
        
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
        mkdir([dirname_DCM 'PEB_firstlevel' filesep 'tempdir_' p.diagnosis_list{k} '_' conditions{c}])
        cd([dirname_DCM 'PEB_firstlevel' filesep 'tempdir_' p.diagnosis_list{k} '_' conditions{c}])
        [PEB,DCM] = spm_dcm_peb(DCM,M,'all');
        DCM = spm_dcm_reduce(DCM,rE,rC);
        cd(this_dir)
        if ~exist([dirname_DCM 'PEB_firstlevel' filesep])
            mkdir([dirname_DCM 'PEB_firstlevel' filesep])
        end
        
        save([dirname_DCM 'PEB_firstlevel' filesep 'First_level_PEB_' p.diagnosis_list{k} '_' conditions{c} '.mat'],'DCM','PEB')
    end
end

    function restore_env(old_path)
        path(old_path);
    end

end