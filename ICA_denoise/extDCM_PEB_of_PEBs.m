function extDCM_PEB_of_PEBs(dirname_DCM,conditions,groups,p,PEB_focus,regions,conductances)

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

assert(length(PEB_focus)==1,'More than one PEB focus specified, which does not work for this script')

for k = groups
    these_labels = regexp(PEB_focus(1),'\d*','Match');
    if isempty(these_labels{1})
        this_tail = PEB_focus{1};
    elseif length(these_labels{1})==2
        this_tail = [PEB_focus{1}(1) '_' regions{str2num(these_labels{1}{1})} '_' conductances{str2num(these_labels{1}{2})}];
    else
        error('Neither zero nor two numbers in the PEB focus - I do not know how to interpret this')
    end
    
    this_PEB = load(['PEB_' this_tail '_' p.diagnosis_list{k} '.mat'],'PEB');
    PEB{k} = this_PEB.PEB;
end

X = ones(length(groups),1); %Group mean
for i = 1:length(groups) %Every pairwise contrast once
    j = i;
    while j<=length(groups)-1
        j = j+1;
        X(i,end+1) = 1;
        X(j,end) = -1;
    end
end

mkdir([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_' this_tail '_' char(join(conditions,'_')) '_Overall'])
cd([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_' this_tail '_' char(join(conditions,'_')) '_Overall'])
[PEB_Overall,DCM_Overall] = spm_dcm_peb(PEB',X,PEB_focus(1));
cd([dirname_DCM 'PEB_secondlevel'])
save(['PEB_' this_tail '_Overall.mat'],'PEB_Overall','DCM_Overall')

cd([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_' this_tail '_' char(join(conditions,'_')) '_Overall'])
% find the posterior p-vals for the parameters in this PEB
BMA_Overall = spm_dcm_peb_bmc(PEB_Overall);
BMA_Overall.Ep = reshape(BMA_Overall.Ep,[],size(X,2));
BMA_Overall.Cp = reshape(BMA_Overall.Cp,[],size(X,2));
BMA_Overall.Pp = reshape(BMA_Overall.Pp,[],size(X,2));
cd([dirname_DCM 'PEB_secondlevel'])
save(['PEB_' this_tail '_Overall.mat'],'BMA_Overall','-append')

    function restore_env(old_path)
        path(old_path);
    end

end

