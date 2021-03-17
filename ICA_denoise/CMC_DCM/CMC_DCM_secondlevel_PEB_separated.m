function CMC_DCM_secondlevel_PEB_separated(dirname_DCM,conditions,groups,p,PEB_focus,regions,model_number,Participant)

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

template_X = ones(length(groups),1); %Group mean
for i = 1:length(groups) %Every pairwise contrast once
    j = i;
    while j<=length(groups)-1
        j = j+1;
        template_X(i,end+1) = 1;
        template_X(j,end) = -1;
    end
end
X = [];

for c = 1:length(conditions)
    for this_group_number = 1:length(groups)
        k = groups(this_group_number);
        this_group = find(p.group==k);
        if isempty(this_group)
            if strcmp(p.diagnosis_list{k},'All_AD')
                this_group = [find(p.group==find(contains(p.diagnosis_list,'pca')));find(p.group==find(contains(p.diagnosis_list,'ADMCI')))];
            elseif strcmp(p.diagnosis_list{k},'All_FTD')
                this_group = [find(p.group==find(contains(p.diagnosis_list,'bvFTD')));find(p.group==find(contains(p.diagnosis_list,'nfvppa')))];
            end
        end
        dcm_files={};
        for subj = 1:length(this_group)
            tmp = dir([dirname_DCM 'mod_' num2str(model_number) '_*' Participant{this_group(subj)}.namepostmerge '*' conditions{c} '*.mat']);
            try
                dcm_files{end+1} = [dirname_DCM tmp.name];
            catch
                error([dirname_DCM 'mod_' num2str(model_number) '_*' Participant{this_group(subj)}.namepostmerge '*' conditions{c} '*.mat not found'])
            end
        end
        
        save(['dcm_files'  p.diagnosis_list{k} '_' conditions{c}],'dcm_files')
        
        X = [X; repmat(template_X(this_group_number,:),length(this_group),1)];
        
        %  Load the DCMs into a cell array and set up the variables for the PEB:
        clear DCM
        for i = 1:length(dcm_files)
            GCM{end+1,1} = load(dcm_files{i});
            GCM{end,1} = GCM{end,1}.DCM;
        end
    end
    M = GCM{1}.M;
    M.X = X;
    [PEB, DCM2] = spm_dcm_peb(GCM,M, PEB_focus);
    
    mkdir([dirname_DCM 'PEB_secondlevel'])
    save([dirname_DCM 'PEB_secondlevel/PEB_' PEB_focus{1} '_' cat(2,p.diagnosis_list{groups}) '.mat'],'PEB','DCM2')
    
    BMA = spm_dcm_peb_bmc(PEB);
    
    spm_dcm_peb_review(BMA,GCM)
    
    BMA.Ep = reshape(BMA.Ep,[],size(X,2));
    BMA.Cp = reshape(BMA.Cp,[],size(X,2));
    BMA.Pp = reshape(BMA.Pp,[],size(X,2));
    mkdir([dirname_DCM 'PEB_secondlevel'])
    save([dirname_DCM 'PEB_secondlevel/PEB_' PEB_focus{1} '_' cat(2,p.diagnosis_list{groups}) '.mat'],'BMA','-append')

end
        