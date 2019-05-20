% List of open inputs
% Head model specification: M/EEG datasets - cfg_files
% Head model specification: Individual structural image - cfg_files
nrun = X; % enter the number of runs here
jobfile = {'/group/language/data/thomascope/MMN/ICA_denoise/scratch_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(2, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Head model specification: M/EEG datasets - cfg_files
    inputs{2, crun} = MATLAB_CODE_TO_FILL_INPUT; % Head model specification: Individual structural image - cfg_files
end
spm('defaults', 'EEG');
spm_jobman('run', jobs, inputs{:});
