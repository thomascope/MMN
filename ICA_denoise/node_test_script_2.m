function node_test_script_2
% 
% spmpath = '/group/language/data/thomascope/MMN/ICA_denoise/Tallie_extDCM/spm12_latestTA/';
% thisspm = which('spm');
% if ~strcmp(thisspm(1:end-5), spmpath)
%     rmpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/'));
%     rmpath(genpath('/group/language/data/thomascope/spm12_fil_r6906/'));
%     addpath(spmpath)
%     spm eeg
% end

addpath('/group/language/data/thomascope/MMN/ICA_denoise/Tallie_extDCM/mfiles_also_needed');

load('node_test_params.mat','p','all_names')
p = p; %Stupid re-naming required for parfor loop
all_names = all_names;

dirname_DCM = '/imaging/mlr/users/tc02/Holly_MMN/extDCMs_hE6/';
filestem = 'b8LFP_s_-100_500_LOR_fmbraedfffM';
conditions = {'STD','DVT','location','intensity','duration','gap','frequency'};
all_combinations = combvec(unique(p.group)',1:length(conditions));
Poolinfo = cbupool(35,'--mem-per-cpu=8G --time=167:00:00 --nodelist=node-i[06-07]'); %Edit this to whichever two nodes-i are idle
parpool(Poolinfo,Poolinfo.NumWorkers,'SpmdEnabled',false);
try
    parfor this_comb = 1:length(all_combinations)
        %for this_comb = 1:length(all_combinations) %falls over in parallel due to tmp.mat and unpredictable cd behaviour - needs fixing for bigger datasets
        k = all_combinations(1,this_comb)
        c = all_combinations(2,this_comb)
        extDCM_nodetest(dirname_DCM,filestem,conditions(c),k,p,all_names)
    end
catch
    rmdir([pwd filesep 'nodetest_folders' filesep 'tempdir_*'])
    
    delete(gcp)
end

    function restore_env(old_path)
        path(old_path);
    end

end