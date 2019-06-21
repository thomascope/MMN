prefix = 'wrmtf_braedfffM*.mat';
TFweightedgrandaveragecomplete = zeros(1,1);
this_output_folder_tail = {};
figure
for todonumber = 1:nsubj
    if iscell(Participant{todonumber}.name)
        this_output_folder_tail{todonumber}  = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.namepostmerge '/'];
    else
        this_output_folder_tail{todonumber}  = [Participant{todonumber}.groupfolder '/' Participant{todonumber}.name '/'];
    end
    fullpath = dir([pathstem this_output_folder_tail{todonumber} prefix]);
    this_temp_data = spm_eeg_load([pathstem this_output_folder_tail{todonumber} fullpath.name]);
    imagesc(flipud(squeeze(mean(this_temp_data(this_temp_data.selectchannels('MEGPLANAR'),:,find(this_temp_data.time>=p.windows(1,1)/1000&this_temp_data.time<=p.windows(1,2)/1000),1),1))))
    pause
end
    

these_diags = unique(all_diags,'stable');    
figure
colormap('jet')
for i = 1:length(these_files)
subplot(5,1,i)
this_temp_data = spm_eeg_load([these_diags{i} '_grandmean.mat']);
%imagesc(flipud(squeeze(mean(this_temp_data(this_temp_data.selectchannels('MEGPLANAR'),:,find(this_temp_data.time>=p.windows(1,1)/1000&this_temp_data.time<=p.windows(1,2)/1000),1),1))))
imagesc(flipud(squeeze(mean(this_temp_data(this_temp_data.selectchannels('MEGPLANAR'),find(this_temp_data.frequencies<=40),find(this_temp_data.time>=p.windows(1,1)/1000&this_temp_data.time<=p.windows(1,2)/1000),1),1)))-flipud(squeeze(mean(this_temp_data(this_temp_data.selectchannels('MEGPLANAR'),find(this_temp_data.frequencies<=40),find(this_temp_data.time>=p.windows(1,1)/1000&this_temp_data.time<=p.windows(1,2)/1000),2),1))),[-0.1 0.1])
end