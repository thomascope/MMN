%-----------------------------------------------------------------------
% Job saved on 03-Mar-2021 21:59:27 by cfg_util (rev $Rev: 6460 $)
% spm SPM - SPM12 (6906)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.dcm.bms.inference.dir = '<UNDEFINED>';
matlabbatch{1}.spm.dcm.bms.inference.sess_dcm{1}.dcmmat = {
                                                           '/imaging/mlr/users/tc02/Holly_MMN/CMC_DCMs/mod_10_b8LFP_s_-100_500_LOR_fmcffbeMbvftd10_meg14_0114_vc4_dcm_STD_DVT.mat'
                                                           '/imaging/mlr/users/tc02/Holly_MMN/CMC_DCMs/mod_10_b8LFP_s_-100_500_LOR_fmcffbeMbvftd11_meg13_0506_cc410179_dcm_STD_DVT.mat'
                                                           '/imaging/mlr/users/tc02/Holly_MMN/CMC_DCMs/mod_10_b8LFP_s_-100_500_LOR_fmcffbeMbvftd12_meg13_0403_cc510639_dcm_STD_DVT.mat'
                                                           '/imaging/mlr/users/tc02/Holly_MMN/CMC_DCMs/mod_10_b8LFP_s_-100_500_LOR_fmcffbeMbvftd13_meg13_0393_cc520065_dcm_STD_DVT.mat'
                                                           '/imaging/mlr/users/tc02/Holly_MMN/CMC_DCMs/mod_10_b8LFP_s_-100_500_LOR_fmcffbeMbvftd14_meg13_0480_cc420261_dcm_STD_DVT.mat'
                                                           };
matlabbatch{1}.spm.dcm.bms.inference.sess_dcm{2}.dcmmat = '<UNDEFINED>';
matlabbatch{1}.spm.dcm.bms.inference.model_sp = {''};
matlabbatch{1}.spm.dcm.bms.inference.load_f = {''};
matlabbatch{1}.spm.dcm.bms.inference.method = 'RFX';
matlabbatch{1}.spm.dcm.bms.inference.family_level.family_file = {''};
matlabbatch{1}.spm.dcm.bms.inference.bma.bma_no = 0;
matlabbatch{1}.spm.dcm.bms.inference.verify_id = 1;



function jobfile = create_CMC_DCM_BMS_Job(output_directory,input_directory,Participant,group)
%A function for creating an appropriate DCM BMS for a given group

% Work out which subjects to inlcude
if strcmp(group,'all')
    include_this = true(1,length(Participant));
else
    for i = 1:length(Participant)
        include_this(i) = strcmp(Participant{i}.diag,group);
    end
end

Participant = Participant(include_this);


jobfile = ['./BMS_DCM_jobfiles/BMS_' group '_' num2str(date) '_job.m'];
fileID = fopen(jobfile,'w');
%Output directory
fprintf(fileID,['matlabbatch{1}.spm.dcm.bms.inference.dir = {''' char(output_directory) '''}\n']);

%Now define the scans
fprintf(fileID,['matlabbatch{1}.spm.dcm.bms.inference.sess_dcm{1}.dcmmat = {\n']);
for i = 1:length(Participant)
    for j = 1:size(dir([input_directory '*' Participant{1}.namepostmerge '*']),1)
        this_file = dir([input_directory 'mod_' num2str(j) '_*' Participant{1}.namepostmerge '*'])
        fprintf(fileID,['' this_file.name '''\n']

stim_duration = '1.0';

for runI = 1:size(starttime,2)
    %First define the scan locations
    fprintf(fileID,['matlabbatch{1}.spm.stats.fmri_spec.sess(' num2str(runI) ').scans = {\n']);
    for this_scan = 1:size(inputs{(2*(runI-1))+2},1)
        fprintf(fileID,['''' char(inputs{(2*(runI-1))+2}(this_scan,:)) '''\n']);
    end
    fprintf(fileID,['};\n']);
    
    %Then define the condition labels
    for this_cond = 1:size(stim_type_labels{runI},1)
        fprintf(fileID,['matlabbatch{1}.spm.stats.fmri_spec.sess(' num2str(runI) ').cond(' num2str(this_cond) ').name = ''' stim_type_labels{runI}{this_cond,1} '_' stim_type_labels{runI}{this_cond,2} '_' stim_type_labels{runI}{this_cond,3} '_' stim_type_labels{runI}{this_cond,4} ''';\n']);
        fprintf(fileID,['matlabbatch{1}.spm.stats.fmri_spec.sess(' num2str(runI) ').cond(' num2str(this_cond) ').onset = [' num2str(starttime{runI}(stimType{runI}==this_cond)) '];\n']);
        fprintf(fileID,['matlabbatch{1}.spm.stats.fmri_spec.sess(' num2str(runI) ').cond(' num2str(this_cond) ').duration = ' stim_duration ';\n']);
        fprintf(fileID,['matlabbatch{1}.spm.stats.fmri_spec.sess(' num2str(runI) ').cond(' num2str(this_cond) ').tmod = 0;\n']);
        fprintf(fileID,['matlabbatch{1}.spm.stats.fmri_spec.sess(' num2str(runI) ').cond(' num2str(this_cond) ').pmod = struct(''name'', {}, ''param'', {}, ''poly'', {});\n']);
        fprintf(fileID,['matlabbatch{1}.spm.stats.fmri_spec.sess(' num2str(runI) ').cond(' num2str(this_cond) ').orth = 1;\n']);
    end
    %Now model the button press
    fprintf(fileID,['matlabbatch{1}.spm.stats.fmri_spec.sess(' num2str(runI) ').cond(' num2str(this_cond+1) ').name = ''Left Button Press'';\n']);
    fprintf(fileID,['matlabbatch{1}.spm.stats.fmri_spec.sess(' num2str(runI) ').cond(' num2str(this_cond+1) ').onset = [' num2str(buttonpresstime{runI}(buttonpressed{runI}==1)) '];\n']);
    fprintf(fileID,['matlabbatch{1}.spm.stats.fmri_spec.sess(' num2str(runI) ').cond(' num2str(this_cond+1) ').duration = 0;\n']);
    fprintf(fileID,['matlabbatch{1}.spm.stats.fmri_spec.sess(' num2str(runI) ').cond(' num2str(this_cond+1) ').tmod = 0;\n']);
    fprintf(fileID,['matlabbatch{1}.spm.stats.fmri_spec.sess(' num2str(runI) ').cond(' num2str(this_cond+1) ').pmod = struct(''name'', {}, ''param'', {}, ''poly'', {});\n']);
    fprintf(fileID,['matlabbatch{1}.spm.stats.fmri_spec.sess(' num2str(runI) ').cond(' num2str(this_cond+1) ').orth = 1;\n']);
    fprintf(fileID,['matlabbatch{1}.spm.stats.fmri_spec.sess(' num2str(runI) ').cond(' num2str(this_cond+2) ').name = ''Right Button Press'';\n']);
    fprintf(fileID,['matlabbatch{1}.spm.stats.fmri_spec.sess(' num2str(runI) ').cond(' num2str(this_cond+2) ').onset = [' num2str(buttonpresstime{runI}(buttonpressed{runI}==4)) '];\n']);
    fprintf(fileID,['matlabbatch{1}.spm.stats.fmri_spec.sess(' num2str(runI) ').cond(' num2str(this_cond+2) ').duration = 0;\n']);
    fprintf(fileID,['matlabbatch{1}.spm.stats.fmri_spec.sess(' num2str(runI) ').cond(' num2str(this_cond+2) ').tmod = 0;\n']);
    fprintf(fileID,['matlabbatch{1}.spm.stats.fmri_spec.sess(' num2str(runI) ').cond(' num2str(this_cond+2) ').pmod = struct(''name'', {}, ''param'', {}, ''poly'', {});\n']);
    fprintf(fileID,['matlabbatch{1}.spm.stats.fmri_spec.sess(' num2str(runI) ').cond(' num2str(this_cond+2) ').orth = 1;\n']);
    
    
    fprintf(fileID,['matlabbatch{1}.spm.stats.fmri_spec.sess(' num2str(runI) ').multi = {''''};\n']);
    fprintf(fileID,['matlabbatch{1}.spm.stats.fmri_spec.sess(' num2str(runI) ').regress = struct(''name'', {}, ''val'', {});\n']);
    fprintf(fileID,['matlabbatch{1}.spm.stats.fmri_spec.sess(' num2str(runI) ').multi_reg = {''' char(inputs{(2*(runI-1))+3}) '''};\n']);
    fprintf(fileID,['matlabbatch{1}.spm.stats.fmri_spec.sess(' num2str(runI) ').hpf = 128;\n']);
    
end

%Postamble and model estimation
fprintf(fileID,['matlabbatch{1}.spm.stats.fmri_spec.fact = struct(''name'', {}, ''levels'', {});\nmatlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];\nmatlabbatch{1}.spm.stats.fmri_spec.volt = 1;\nmatlabbatch{1}.spm.stats.fmri_spec.global = ''None'';\nmatlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.4;\nmatlabbatch{1}.spm.stats.fmri_spec.mask = {''' char(inputs{end}) '''};\nmatlabbatch{1}.spm.stats.fmri_spec.cvi = ''AR(1)'';\nmatlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep(''fMRI model specification: SPM.mat File'', substruct(''.'',''val'', ''{}'',{1}, ''.'',''val'', ''{}'',{1}, ''.'',''val'', ''{}'',{1}), substruct(''.'',''spmmat''));\nmatlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;\nmatlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;\n']);

% Now create contrasts
fprintf(fileID,['matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep(''Model estimation: SPM.mat File'', substruct(''.'',''val'', ''{}'',{2}, ''.'',''val'', ''{}'',{1}, ''.'',''val'', ''{}'',{1}), substruct(''.'',''spmmat''));\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = ''Photos > Line'';\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = kron([1,-1],ones(1,30)) ;\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = ''bothsc'';\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = ''Photos < Line'';\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = kron([-1,1],ones(1,30)) ;\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = ''bothsc'';\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = ''Left > Right'';\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = repmat(kron([1,-1],ones(1,5)),1,6) ;\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = ''bothsc'';\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = ''Left < Right'';\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = repmat(kron([-1,1],ones(1,5)),1,6) ;\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = ''bothsc'';\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = ''Common > Rare'';\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = repmat(kron([1,0,-1],ones(1,10)),1,2) ;\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = ''bothsc'';\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = ''Common < Rare'';\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = repmat(kron([-1,0,1],ones(1,10)),1,2) ;\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = ''bothsc'';\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = ''Left Button'';\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [zeros(1,60),1] ;\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = ''bothsc'';\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = ''Right Button'';\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [zeros(1,61),1] ;\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = ''bothsc'';\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = ''All Pictures'';\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = repmat(kron([1,1,1],ones(1,10)),1,2) ;\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = ''bothsc'';\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{10}.tcon.name = ''Negative All Pictures'';\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{10}.tcon.weights = repmat(kron([-1,-1,-1],ones(1,10)),1,2) ;\n']);
fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{10}.tcon.sessrep = ''bothsc'';\n']);

%Now create condition contrasts for later RSA
for this_cond = 1:size(stim_type_labels{runI},1)
    fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{' num2str(10+this_cond) '}.tcon.name= ''Condition ' num2str(this_cond) ''';\n']);
    fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{' num2str(10+this_cond) '}.tcon.weights = [' num2str([zeros(1,this_cond-1) 1 zeros(1,size(stim_type_labels{runI},1)-this_cond)]) '];\n']);
    fprintf(fileID,['matlabbatch{3}.spm.stats.con.consess{' num2str(10+this_cond) '}.tcon.sessrep = ''replsc'';\n']);
end

fprintf(fileID,['matlabbatch{3}.spm.stats.con.delete = 0;\n']);

