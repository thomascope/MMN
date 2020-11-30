%% Pre-processing - merge any files recorded as separate runs
parfor cnt = 1:length(These_Participants)
    Preprocessing_mainfunction('merge','secondfilter',p,[pathstem These_Participants{cnt}.groupfolder '/'], [], These_Participants{cnt}.namepostmerge,cnt,[],[],[],[], badeeg)
end
%% Pre-processing - sort trial order
parfor cnt = 1:length(These_Participants)
    Preprocessing_mainfunction('sort','merge',p,[pathstem These_Participants{cnt}.groupfolder '/'], [], These_Participants{cnt}.namepostmerge,cnt,[],[],[],[], badeeg)
end
%% Now specify forward model
re_forwardmodelcomplete = zeros(1,length(These_Participants));
parfor todonumber = 1:length(These_Participants)
    if re_forwardmodelcomplete(todonumber)~=1
        
        try
            These_Participants{todonumber}.name = These_Participants{todonumber}.namepostmerge;
        end
        megpaths = {[pathstem These_Participants{todonumber}.groupfolder '/' These_Participants{todonumber}.name '/fmcffbeM' These_Participants{todonumber}.name '.mat'],
            [pathstem These_Participants{todonumber}.groupfolder '/' These_Participants{todonumber}.name '/cffbeM' These_Participants{todonumber}.name '.mat']
            };
        if ~exist(megpaths{1},'file')
            for i = 1:length(megpaths)
                megpaths{i} = ls([megpaths{i}(1:end-4) '*.mat'])
            end
        end
        mripath = [mridirectory These_Participants{todonumber}.groupfolder '/' These_Participants{todonumber}.name '/' These_Participants{todonumber}.MRI '.nii'];
        if ~exist(mripath,'file') && strcmp(These_Participants{todonumber}.MRI,'single_subj_T1')
            mripath = ['/group/language/data/thomascope/spm12_fil_r6906/canonical/single_subj_T1.nii'];
        end
        newmripath = [pathstem These_Participants{todonumber}.groupfolder '/' These_Participants{todonumber}.name '/MRI/'];
        if ~exist(newmripath,'dir')
            mkdir(newmripath)
        end
        try
            copyfile(mripath, newmripath)
        catch
            mripath = ['/group/language/data/thomascope/spm12_fil_r6906/canonical/single_subj_T1.nii'];
            warning(['Scan missing for subject ' num2str(todonumber) ' using template instead.']);
            copyfile(mripath, newmripath);
            These_Participants{todonumber}.MRI = 'single_subj_T1';
        end
        newmripath = [pathstem These_Participants{todonumber}.groupfolder '/' These_Participants{todonumber}.name '/MRI/' These_Participants{todonumber}.MRI '.nii'];
        try
            forward_model_this_subj(megpaths,newmripath, p)
            re_forwardmodelcomplete(todonumber) = 1;
            fprintf('\n\nForward modelling complete for subject number %d,\n\n',todonumber);
        catch
            re_forwardmodelcomplete(todonumber) = 0;
            fprintf('\n\nForward modelling failed for subject number %d\n\n',todonumber);
        end
    end
end
