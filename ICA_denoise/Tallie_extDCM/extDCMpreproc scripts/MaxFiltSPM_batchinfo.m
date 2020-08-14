%% Set up global variables

% clear all

% make sure EEG modality of SPM software is selected
spm('defaults','EEG');

% define paths
rawpathstem = '/megdata/cbu/ftdrug/';
pathstem = '/imaging/na01/MEGmf_mem/';

% define conditions
conditions = {'Dev' 'rep1' 'rep2' 'rep3' 'rep4' 'rep5' 'rep6' 'rep7' 'rep8' 'rep9'};

contrast_labels = {};
contrast_weights = [];    

subjects   = {'' ''}; % subject folder on the CBUs /megdata
dates      = {'' ''}; % date folder inside the subject folder
blocksin   = {{''} {''}}; % orig file name (without the raw suffix)
blocksout  = {{''} {''}}; % what you want it to be called

load('/imaging/na01/misc/extDCMstuff/bad_electrodes_mem.mat')
badchannels = {};
badeeg = {};
for k = 1:length(subjects)
    for k2 = 2:3
        i = find(strcmp(subjects{k},te(:,k2)));
        if ~isempty(i)
            if ~isempty(te{i,k2+2})
                badchannels{k,1} = te{i,k2+2};
            end
            if ~isempty(te{i,k2+4})
                badeeg{k,1} = te{i,k2+4};
            end
        end
    end
end
group               = 1 + cellfun(@(x) isempty(regexpi(x,'p','once')),blocksin);

for i = 1:length(subjects) % cnt
    if strcmp(class(subjects{i}),'char')
        subjects1{i}.newname = subjects{i};
        subjects1{i}.oldname = subjects{i};
    else
        subjects1{i}.newname = subjects{i}.newname;
        subjects1{i}.oldname = subjects{i}.oldname;
    end
end
subjects = subjects1;   

