function extDCM_PEB_of_PEBs(dirname_DCM,conditions,groups,p,PEB_focus,regions,conductances,combinegroups)

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

these_labels = regexp(PEB_focus(1),'\d*','Match');
if isempty(these_labels{1})
    this_tail = PEB_focus{1};
elseif length(these_labels{1})==2
    this_tail = [PEB_focus{1}(1) '_' regions{str2num(these_labels{1}{1})} '_' conductances{str2num(these_labels{1}{2})}];
else
    error('Neither zero nor two numbers in the PEB focus - I do not know how to interpret this')
end

if combinegroups == 0
  for k = groups
   
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
elseif combinegroups == 1
    PEB = {};
    for k = [1, max(unique(p.group))+1:max(unique(p.group))+2] % Three groups - controls, all ftd, all AD
        this_PEB = load(['PEB_' this_tail '_' p.diagnosis_list{k} '.mat'],'PEB');
        PEB{end+1} = this_PEB.PEB;
    end
    X = ones(3,1); %Group mean
    for i = 1:3 %Every pairwise contrast once
        j = i;
        while j<=3-1
            j = j+1;
            X(i,end+1) = 1;
            X(j,end) = -1;
        end
    end
    
end

if combinegroups == 0
mkdir([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_' this_tail '_' char(join(conditions,'_')) '_Overall'])
cd([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_' this_tail '_' char(join(conditions,'_')) '_Overall'])
elseif combinegroups == 1
    mkdir([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_' this_tail '_' char(join(conditions,'_')) '_Overall_combined'])
cd([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_' this_tail '_' char(join(conditions,'_')) '_Overall_combined'])
end
[PEB_Overall,DCM_Overall] = spm_dcm_peb(PEB',X,PEB_focus(1));
cd([dirname_DCM 'PEB_secondlevel'])
if combinegroups == 0
    save(['PEB_' this_tail '_Overall.mat'],'PEB_Overall','DCM_Overall')
    cd([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_' this_tail '_' char(join(conditions,'_')) '_Overall'])
elseif combinegroups == 1
    save(['PEB_' this_tail '_Overall_combined.mat'],'PEB_Overall','DCM_Overall')
    cd([dirname_DCM 'PEB_secondlevel' filesep 'tempdir_' this_tail '_' char(join(conditions,'_')) '_Overall_combined'])
end

% find the posterior p-vals for the parameters in this PEB
BMA_Overall = spm_dcm_peb_bmc(PEB_Overall);
BMA_Overall.Ep = reshape(BMA_Overall.Ep,[],size(X,2));
BMA_Overall.Cp = reshape(BMA_Overall.Cp,[],size(X,2));
BMA_Overall.Pp = reshape(BMA_Overall.Pp,[],size(X,2));
cd([dirname_DCM 'PEB_secondlevel'])
if combinegroups == 0
    save(['PEB_' this_tail '_Overall.mat'],'BMA_Overall','-append')
elseif combinegroups == 1
    save(['PEB_' this_tail '_Overall_combined.mat'],'BMA_Overall','-append')
end

    function restore_env(old_path)
        path(old_path);
    end

end

