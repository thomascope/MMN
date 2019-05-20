function forward_model_this_subj(megpath, mripath)
%forward models a cellstring of megpaths against a cellstring of mripaths

% List of open inputs
% Head model specification: M/EEG datasets - cfg_files
% Head model specification: Individual structural image - cfg_files
nrun = 1; % enter the number of runs here
%jobfile = {'/group/language/data/thomascope/vespa/SPM12version/Standalone preprocessing pipeline/tc_source/batch_forwardmodel_job_noheadpoints.m'};
jobfile = {[pwd '/batch_forward_auto.m']};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(5, nrun);
mrilist = mripath;
meglist = megpath;


load('/imaging/tc02/vespa/preprocess/controls_mni');
load('/imaging/tc02/vespa/preprocess/patients_mni');
%allmni = [corrected_controlsmni, corrected_patientsmni]';
allmni = [controlsmni, patientsmni]';
                                             
% if length(meglist)~=nrun || length(mrilist)~=nrun || length(datadir_early)~=nrun || length(datadir_late)~=nrun
%     error('The number of runs specified is not the same as the length of one of the inputs. Please double check this.')
% end
for crun = 1:nrun
    inputs{1, crun} = cellstr(meglist{crun}); % Head model specification: M/EEG datasets - cfg_files
    inputs{2, crun} = cellstr(mrilist{crun}); % Head model specification: Individual structural image - cfg_files
end

%% Now compute the forwards model
forwardmodelworkedcorrectly = zeros(1,nrun);
jobs = repmat(jobfile, 1, 1);

parfor crun = 1:nrun
    spm('defaults', 'EEG');
    spm_jobman('initcfg')
    try
        spm_jobman('run', jobs, inputs{:,crun});
        forwardmodelworkedcorrectly(crun) = 1;
    catch
        forwardmodelworkedcorrectly(crun) = 0;
    end
end