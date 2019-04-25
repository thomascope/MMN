function [D,new_filepath] = DoInvTA_TEC(D,folder,outfolder, varargin)

[f1,f2,f3] = fileparts(D);
if isempty(varargin), varargin{1} = 1; end
new_filepath = sprintf('%s/%s/LOR_%s%s',outfolder,folder, f2, f3);

S.D = D;
S.outfile = new_filepath;
if ~exist(S.outfile,'file')
spm_eeg_copy(S);
end

if isempty(varargin), varargin{1} = 1; end
Meth = {'IID','LOR'};

m{1}.spm.meeg.source.invert.D = {new_filepath};
m{1}.spm.meeg.source.invert.val = varargin{1};
m{1}.spm.meeg.source.invert.whatconditions.all = 1;
m{1}.spm.meeg.source.invert.isstandard.custom.invtype = Meth{varargin{1}};
m{1}.spm.meeg.source.invert.isstandard.custom.woi = [-Inf Inf];
m{1}.spm.meeg.source.invert.isstandard.custom.foi = [0 256];
m{1}.spm.meeg.source.invert.isstandard.custom.hanning = 1;
m{1}.spm.meeg.source.invert.isstandard.custom.priors.priorsmask = {''};
m{1}.spm.meeg.source.invert.isstandard.custom.priors.space = 1;
m{1}.spm.meeg.source.invert.isstandard.custom.restrict.locs = zeros(0, 3);
m{1}.spm.meeg.source.invert.isstandard.custom.restrict.radius = 32;
m{1}.spm.meeg.source.invert.isstandard.custom.restrict.mask = {''};
m{1}.spm.meeg.source.invert.modality = {'MEGPLANAR','MEG'}; %comment out for multimodal

spm('defaults','EEG');
spm_jobman('run',m);

% [f1,f2,f3] = fileparts(D);
% load(D)
% D.fname = ['I' f2 f3];
% D.data.fname = [f1 filesep 'I' f2 '.dat'];
% save([f1 filesep 'I' f2 f3],'D')
% delete([f1 filesep f2 f3])
% movefile([f1 filesep f2 '.dat'],[f1 filesep 'I' f2 '.dat'])


end

