function D = PP_Coreg(D,struct)


od = cd(fileparts(D));
spm_jobman('initcfg'); 

% specify:    
m{1}.spm.meeg.source.headmodel.D = {D};
m{1}.spm.meeg.source.headmodel.val = 1;
m{1}.spm.meeg.source.headmodel.comment = '';
if strcmp('',struct)
    m{1}.spm.meeg.source.headmodel.meshing.meshes.template = 1; warning('!!! --- USING TEMPLATE STRUCTURAL --- !!!')
else m{1}.spm.meeg.source.headmodel.meshing.meshes.mri = {struct};
end
m{1}.spm.meeg.source.headmodel.meshing.meshres = 2;
m{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'Nasion';
m{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.select = 'nas';
m{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'LPA';
m{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.select = 'lpa';
m{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'RPA';
m{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.select = 'rpa';
m{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 1;
m{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
m{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';

spm('defaults','EEG');
spm_jobman('run',m);

% [f1,f2,f3] = fileparts(D);
% load(D)
% D.fname = ['C' f2 f3];
% D.data.fname = [f1 filesep 'C' f2 '.dat'];
% save([f1 filesep 'C' f2 f3],'D')
% delete([f1 filesep f2 f3])
% movefile([f1 filesep f2 '.dat'],[f1 filesep 'C' f2 '.dat'])

cd(od)


end

