% List of open inputs
% Model Inference: Models - cfg_files
nrun = X; % enter the number of runs here
jobfile = {'/group/language/data/thomascope/MMN/ICA_denoise/CMC_DCM/BMS_CMC_DCM_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(1, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Model Inference: Models - cfg_files
end
spm('defaults', 'EEG');
spm_jobman('run', jobs, inputs{:});
