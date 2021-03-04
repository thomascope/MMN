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

if ~exist('./BMS_DCM_jobfiles')
    mkdir('./BMS_DCM_jobfiles')
end
if ~exist(output_directory)
    mkdir(output_directory)
end
jobfile = ['./BMS_DCM_jobfiles/BMS_' group '_' num2str(date) '_job.m'];
fileID = fopen(jobfile,'w');
%Output directory
fprintf(fileID,['matlabbatch{1}.spm.dcm.bms.inference.dir = {''' char(output_directory) '''}\n']);

%Now define the scans
for i = 1:length(Participant)
    fprintf(fileID,['matlabbatch{1}.spm.dcm.bms.inference.sess_dcm{' num2str(i) '}.dcmmat = {\n']);
    for j = 1:size(dir([input_directory '*' Participant{i}.namepostmerge '*']),1)
        this_file = dir([input_directory 'mod_' num2str(j) '_*' Participant{i}.namepostmerge '*']);
        fprintf(fileID,['''' input_directory this_file.name '''\n']);
    end
    fprintf(fileID,['};\n']);
    disp(['Done Participant ' num2str(i) ' of ' num2str(length(Participant))])
end

fprintf(fileID,['matlabbatch{1}.spm.dcm.bms.inference.model_sp = {''''};\n']);
fprintf(fileID,['matlabbatch{1}.spm.dcm.bms.inference.load_f = {''''};\n']);
fprintf(fileID,['matlabbatch{1}.spm.dcm.bms.inference.method = ''RFX'';\n']);
fprintf(fileID,['matlabbatch{1}.spm.dcm.bms.inference.family_level.family_file = {''''};\n']);
fprintf(fileID,['matlabbatch{1}.spm.dcm.bms.inference.bma.bma_no = 0;\n']);
fprintf(fileID,['matlabbatch{1}.spm.dcm.bms.inference.verify_id = 1;\n']);


