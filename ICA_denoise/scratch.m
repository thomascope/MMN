for c=1:length(conditions)
    if strcmp(modality{m},'EEG')
        if rejecteeg{s} == 1
            %files{p.group(s)}{s}{c} = [];
        else
            files{p.group(s)}{s}{c} = strjoin([pathstem subjects{s} '/' modality{m} filetype '/' imagetype 'condition_' conditions{c} '.nii'],'');
            if exist(files{p.group(s)}{s}{c},'file')
            else
                files{p.group(s)}{s}{c} = strjoin([pathstem subjects{s} '/' modality{m} filetypesplit '/' imagetype 'condition_' conditions{c} '.nii'],'');
            end
        end
        
    else
        files{p.group(s)}{s}{c} = strjoin([pathstem subjects{s} '/' modality{m} filetype '/' imagetype 'condition_' conditions{c} '.nii'],'');
        if exist(files{p.group(s)}{s}{c},'file')
        else
            files{p.group(s)}{s}{c} = strjoin([pathstem subjects{s} '/' modality{m} filetypesplit '/' imagetype 'condition_' conditions{c} '.nii'],'');
        end
    end
end